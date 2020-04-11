#include <algorithm>
#include <chrono>
#include <string>
#include <vector>

// util
#include "dice_roll.hpp"
#include "pseudo_random.hpp"
#include "stats.hpp"

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
using login_atom = caf::atom_constant<caf::atom("login")>;
using finished_atom = caf::atom_constant<caf::atom("finished")>;
using poke_atom = caf::atom_constant<caf::atom("poke")>;
using disconnect_atom = caf::atom_constant<caf::atom("disconnect")>;
using confirm_atom = caf::atom_constant<caf::atom("confirm")>;
using print_atom = caf::atom_constant<caf::atom("print")>;
using collect_atom = caf::atom_constant<caf::atom("collect")>;
using apply_atom = caf::atom_constant<caf::atom("apply")>;
using complete_atom = caf::atom_constant<caf::atom("complete")>;
using append_atom = caf::atom_constant<caf::atom("append")>;

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
                   const caf::actor initiator) {
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
  caf::actor directory;
  dice_roll dice;
  pseudo_random rand;
};

caf::behavior client(caf::stateful_actor<client_state>* self,
                     const std::uint64_t id, const caf::actor directory,
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
          break;
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

struct directory_state {
  client_map clients;
  pseudo_random random;
  std::uint32_t befriend;
  bool is_poker = false;
  caf::actor poker;
};

caf::behavior directory(caf::stateful_actor<directory_state>* self,
                        std::uint64_t seed, std::uint32_t befriend) {
  auto& s = self->state;
  s.random = pseudo_random(seed);
  s.befriend = befriend;
  return {
    [=](login_atom, std::uint64_t id) {
      auto& s = self->state;
      auto new_client = self->spawn(client, id, self, s.random.next_int());
      s.clients.emplace(id, new_client);
      for (auto& client : s.clients) {
        if (s.random.next_int(100) < s.befriend) {
          self->send(client.second, befriend_atom::value, new_client);
          self->send(new_client, befriend_atom::value, client.second);
        }
      }
    },
    [=](login_atom, std::uint64_t id) {
      self->send(self->state.clients.at(id), logout_atom::value);
    },
    [=](left_atom, std::uint64_t id) {
      auto& s = self->state;
      s.clients.erase(id);
      if (s.clients.empty() && s.is_poker) {
        self->send(s.poker, finished_atom::value);
      }
    },
    [=](poke_atom, behavior_factory& behavior, const caf::actor& accumulator) {
      for (auto& client : self->state.clients) {
        self->send(client.second, act_atom::value, behavior, accumulator);
      }
    },
    [=](disconnect_atom, caf::actor& poker) {
      auto& s = self->state;
      s.is_poker = true;
      s.poker = poker;

      for (auto& client : s.clients) {
        self->send(client.second, logout_atom::value);
      }
    },
  };
}

using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;
struct accumulator_state {
  caf::actor poker;
  time_point start;
  time_point end;
  /// time in milliseconds
  double duration;
  size_t expected;
  bool did_stop = false;
};

caf::behavior accumulator(caf::stateful_actor<accumulator_state>* self,
                          caf::actor poker, size_t expected) {
  auto& s = self->state;
  s.poker = poker;
  s.start = std::chrono::high_resolution_clock::now();
  s.expected = expected;
  return {
    [=](bump_atom, const size_t expected) {
      auto& s = self->state;
      s.expected = (s.expected + expected) - 1;
    },
    [=](stop_atom) {
      auto& s = self->state;
      s.expected--;
      if (s.expected == 1) {
        s.end = std::chrono::high_resolution_clock::now();
        s.duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                       s.end - s.start)
                       .count();
        s.did_stop = true;

        self->send(s.poker, confirm_atom::value);
      }
    },
    [=](print_atom, const caf::actor& poker, size_t i, size_t j) {
      self->send(poker, collect_atom::value, i, j, self->state.duration);
    },
  };
}

struct poker_state {
  std::uint64_t clients;
  std::size_t logouts;
  std::size_t confirmations;
  std::uint64_t turns;
  std::size_t iteration;
  std::vector<caf::actor> directories;
  std::vector<caf::actor> runtimes;
  std::size_t accumulations;
  std::vector<std::vector<double>> finals;
  behavior_factory factory;
  caf::actor bench;
  bool last;
  std::vector<double> turn_series;
};

caf::behavior poker(caf::stateful_actor<poker_state>* self,
                    std::uint64_t clients, std::uint64_t turns,
                    std::vector<caf::actor>& directories,
                    behavior_factory& factory) {
  auto& s = self->state;
  s.clients = clients;
  s.turns = turns;
  s.directories = directories;
  s.factory = factory;
  return {
    [=](apply_atom, caf::actor& bench, bool last) {
      auto& s = self->state;
      s.confirmations = s.turns;
      s.logouts = s.directories.size();
      s.bench = bench;
      s.last = last;
      s.accumulations = 0;

      std::uint64_t turns = s.turns;
      std::size_t index = 0;
      std::vector<double> values;
      values.reserve(s.turns);
      s.finals.push_back(values);

      for (size_t client = 0; client < s.clients; ++client) {
        index = client % s.directories.size();
        self->send(s.directories.at(index), login_atom::value, client);
      }
      // feetback loop?
      for (; turns >= 0; --turns) {
        auto accu = self->spawn(accumulator, self, s.clients);
        for (auto& directory : s.directories) {
          self->send(directory, poke_atom::value, s.factory, accu);
        }

        s.runtimes.push_back(accu);
      }
    },
    [=](confirm_atom) {
      auto& s = self->state;
      --s.confirmations;
      if (s.confirmations == 1) {
        for (auto& d : s.directories) {
          self->send(d, disconnect_atom::value, self);
        }
      }
    },
    [=](finished_atom) {
      auto& s = self->state;
      --s.logouts;
      if (s.logouts == 1) {
        std::size_t turn = 0;

        for (auto& accumulator : s.runtimes) {
          ++s.accumulations;
          self->send(accumulator, print_atom::value, self, s.iteration, turn);
          ++turn;
        }

        s.runtimes.clear();
      }
    },
    [=](collect_atom, std::size_t i, std::size_t j, double duration) {
      auto& s = self->state;
      s.finals[i][j] = duration;
      s.turn_series.push_back(duration);

      --s.accumulations;
      if (s.accumulations == 1) {
        ++s.iteration;
        self->send(s.bench, complete_atom::value);

        if (s.last) {
          sample_stats stats(s.turn_series);
          std::vector<std::vector<double>> turns;
          std::vector<double> qos;

          // TODO ask pony about line 381 to 391

          for (std::size_t l = 0; l < s.finals.size(); ++l) {
            qos.push_back(sample_stats(s.finals.back()).stddev());
            s.finals.pop_back();
          }

          std::stringstream title_text;
          title_text << std::string(31, ' ') << std::string(12, ' ') << "j-mean"
                     << std::string(10, ' ') << "j-median"
                     << std::string(11, ' ') << "j-error"
                     << std::string(10, ' ') << "j-stddev"
                     << std::string(14, ' ') << "quality of service"
                     << std::endl;

          self->send(s.bench, append_atom::value, title_text.str());

          std::stringstream result_text;
          result_text
            << "Turns" << std::string(27, ' ') << std::setw(18)
            << std::setprecision(std::numeric_limits<double>::digits10)
            << stats.mean() << " " << std::setw(18)
            << std::setprecision(std::numeric_limits<double>::digits10)
            << stats.median() << " " << std::setw(18)
            << std::setprecision(std::numeric_limits<double>::digits10)
            << stats.err() << " " << std::setw(18)
            << std::setprecision(std::numeric_limits<double>::digits10)
            << stats.stddev() << " " << std::setw(18)
            << std::setprecision(std::numeric_limits<double>::digits10)
            << sample_stats(qos).median() << std::endl;

          self->send(s.bench, append_atom::value, result_text.str());
        }
      }
    },
  };
}

struct chatapp_state {
  std::uint64_t clients;
  std::uint64_t turns;
  std::vector<caf::actor> directories;
  behavior_factory factory;
  caf::actor poker;
};

caf::behavior chatapp(caf::stateful_actor<chatapp_state>* self) {
  auto& s = self->state;
  s.clients = 1024;
  s.turns = 20;

  std::size_t directories = 8;
  std::uint64_t compute = 50;
  std::uint64_t post = 80;
  std::uint64_t leave = 25;
  std::uint64_t invite = 25;
  std::uint64_t befriend = 10;
  pseudo_random rand(42);

  s.factory = behavior_factory(compute, post, leave, invite);

  for (std::size_t i = 0; i < directories; ++i) {
    s.directories.emplace_back(
      self->spawn(directory, rand.next_int(), befriend));
  }

  s.poker = self->spawn(poker, s.clients, s.turns, s.directories, s.factory);

  return {[=](apply_atom, caf::actor& async_benchmark_completion, bool last) {
    self->send(self->state.poker, async_benchmark_completion, last);
  }};
}
