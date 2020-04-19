#pragma once

#include "pseudo_random.hpp"

class dice_roll {
public:
  dice_roll()
    : random_(42){
      // nop
    };

  explicit dice_roll(pseudo_random rand) : random_(rand) {
    // nop
  }

  bool apply(std::uint64_t probability) {
    return random_.next_int(100) < static_cast<std::uint32_t>(probability);
  }

private:
  pseudo_random random_;
};
