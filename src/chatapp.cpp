#include <algorithm>
#include <chrono>
#include <vector>

// util
#include "dice_roll.hpp"

// caf
#include "caf/all.hpp"

// types
using client_map = std::unordered_map<std::uint64_t, caf::actor>;
using friend_set = std::unordered_set<caf::actor>;
using chat_set = std::unordered_set<caf::actor>;
using client_set = std::unordered_set<caf::actor>;

/// clients actions
enum action { post, leave, invite, compute, none };

// messages
using post_atom = caf::atom_constant<caf::atom("post")>;
using forward_atom = caf::atom_constant<caf::atom("forward")>;
using bump_atom = caf::atom_constant<caf::atom("bump")>;
using stop_atom = caf::atom_constant<caf::atom("stop")>;
using join_atom = caf::atom_constant<caf::atom("join")>;
using leave_atom = caf::atom_constant<caf::atom("leave")>;
using left_atom = caf::atom_constant<caf::atom("left")>;
using befriend_atom = caf::atom_constant<caf::atom("befriend")>;
using logout_atom = caf::atom_constant<caf::atom("logout")>;
using invite_atom = caf::atom_constant<caf::atom("invite")>;
using act_atom = caf::atom_constant<caf::atom("act")>;

/// simulates extern client events for each turn
struct behavior_factory {
  behavior_factory() = default;

  behavior_factory(std::uint64_t compute, std::uint64_t post,
                   std::uint64_t leave, std::uint64_t invite)
    : compute(compute), post(post), leave(leave), invite(invite) {
    // nop
  }

  behavior_factory(const behavior_factory& f) = default;

  action apply(dice_roll dice) {
    auto next_action = action::none;

    if (dice.apply(compute)) {
      next_action = action::compute;
    } else if (dice.apply(post)) {
      next_action = action::post;
    } else if (dice.apply(leave)) {
      next_action = action::leave;
    } else if (dice.apply(invite)) {
      next_action = action::invite;
    }

    return next_action;
  }

  std::uint64_t compute;
  std::uint64_t post;
  std::uint64_t leave;
  std::uint64_t invite;
};

template <class Inspector>
typename Inspector::result_type inspect(Inspector& f, behavior_factory& x) {
  return f(caf::meta::type_name("behavior_factory"), x.compute, x.post, x.leave,
           x.invite);
}

/// for clients compute turn
std::uint64_t fibonacci(std::uint8_t x) {
  if (x == 0 || x == 1) {
    return x;
  } else {
    auto j = x / 2;
    auto fib_j = fibonacci(j);
    auto fib_i = fibonacci(j - 1);

    if (x % 2 == 0) {
      return fib_j * (fib_j + (fib_i * 2));
    } else if (x % 4 == 1) {
      return (((fib_j * 2) + fib_i) * ((fib_j * 2) - fib_i)) + 2;
    } else {
      return (((fib_j * 2) + fib_i) * ((fib_j * 2) - fib_i)) - 2;
    }
  }
}

struct chat_state {
  client_set members;
  std::vector<std::vector<std::uint8_t>> buffer;
};

caf::behavior chat(caf::stateful_actor<chat_state>* self,
                   const caf::actor& initiator) {
  self->state.members.insert(initiator);
  return {
    [=](post_atom, std::vector<std::uint8_t>& payload,
        const caf::actor& accumulator) {
      auto& s = self->state;
#ifndef BENCH_NO_BUFFERED_CHATS
      s.buffer.push_back(payload);
#endif
      if (s.members.empty()) {
        self->send(accumulator, stop_atom::value);
      } else {
        self->send(accumulator, bump_atom::value, s.members.size());
        for (auto& member : s.members) {
          self->send(member, forward_atom::value, self, payload, accumulator);
        }
      }
    },
    [=](join_atom, const caf::actor& client, const caf::actor& accumulator) {
      auto& s = self->state;
      s.members.insert(client);
#ifdef BENCH_NO_BUFFERED_CHATS
      self->send(accumulator, stop_atom::value);
#else
      if (s.buffer.empty()) {
        self->send(accumulator, stop_atom::value);
      } else {
        self->send(accumulator, bump_atom::value, s.buffer.size());
        for (auto& message : s.buffer) {
          self->send(client, forward_atom::value, message, accumulator);
        }
      }
#endif
    },
    [=](leave_atom, const caf::actor& client, const bool did_logout,
        const caf::actor& accumulator) {
      self->state.members.erase(client);
      self->send(client, left_atom::value, self, did_logout, accumulator);
    },
    [=](leave_atom, const caf::actor& client, const bool did_logout) {
      self->state.members.erase(client);
      self->send(client, left_atom::value, self, did_logout);
    },
  };
}

struct client_state {
  std::uint64_t id;
  friend_set friends;
  chat_set chats;
  caf::actor& directory;
  dice_roll dice;
  pseudo_random rand;
};

caf::behavior client(caf::stateful_actor<client_state>* self,
                     const std::uint64_t id, const caf::actor& directory,
                     std::uint64_t seed) {
  auto& s = self->state;
  s.id = id;
  s.directory = directory;
  s.dice = dice_roll(seed);
  s.rand = pseudo_random(seed);
  return {
    [=](befriend_atom, const caf::actor& client) {
      self->state.friends.insert(client);
    },
    [=](logout_atom) {
      auto& s = self->state;
      if (s.chats.empty()) {
        self->send(s.directory, left_atom::value, s.id);
      } else {
        for (auto& chat : s.chats) {
          self->send(chat, leave_atom::value, self, true);
        }
      }
    },
    [=](left_atom, const caf::actor& chat, const bool did_logout,
        const caf::actor& accumulator) {
      auto& s = self->state;
      s.chats.erase(chat);
      self->send(accumulator, stop_atom::value);
    },
    [=](left_atom, const caf::actor& chat, const bool did_logout,
        const caf::actor& accumulator) {
      auto& s = self->state;
      s.chats.erase(chat);
      if (did_logout && s.chats.empty()) {
        self->send(s.directory, left_atom::value, s.id);
      }
    },
    [=](invite_atom, const caf::actor& chat, const caf::actor& accumulator) {
      self->state.chats.insert(chat);
      self->send(chat, join_atom::value, self, accumulator);
    },
    [=](forward_atom, const caf::actor& chat,
        std::vector<std::uint8_t>& payload, const caf::actor& accumulator) {
      self->send(accumulator, stop_atom::value);
    },
    [=](act_atom, behavior_factory& behavior, const caf::actor& accumulator) {
      auto& s = self->state;
      caf::actor next_chat;
      if (s.chats.empty()) {
        next_chat = self->spawn(chat, self);
      } else {
        size_t index = s.rand.next_int(s.chats.size());
        next_chat = *std::next(s.chats.begin(), index);
      }

      switch (behavior.apply(s.dice)) {
        case action::post:
          self->send(next_chat, post_atom::value, std::vector<std::uint8_t>{},
                     accumulator);
          break;
        case action::leave:
          self->send(next_chat, leave_atom::value, false, accumulator);
          break;
        case action::compute:
          fibonacci(35);
          self->send(accumulator, stop_atom::value);
        case action::none:
          self->send(accumulator, stop_atom::value);
          break;
        case action::invite:
          std::vector<caf::actor> f;
          for (auto& i : s.friends) {
            f.push_back(i);
          }
          auto rng = std::default_random_engine{s.rand.next_int()};
          std::shuffle(f.begin(), f.end(), rng);
          f.insert(f.begin(), self);

          size_t invitations = s.rand.next_int(s.friends.size());
          if (invitations == 0) {
            self->send(accumulator, stop_atom::value);
          } else {
            self->send(accumulator, bump_atom::value, invitations);
            for (size_t i = 0; i < invitations; ++i) {
              self->send(f[i], invite_atom::value, next_chat, accumulator);
            }
          }
          break;
      }
    },
  };
}
