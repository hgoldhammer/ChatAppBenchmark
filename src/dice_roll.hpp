#pragma once

#import "pseudo_random.hpp"

class dice_roll {
public:
  dice_roll(std::uint64_t seed) : random_(seed) {
    //nop
  }

  bool apply(std::uint64_t probability) {
    return random_.next_int(100) < static_cast<std::uint32_t>(probability);
  }

private:
  pseudo_random random_;
};
