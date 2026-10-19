#include "core/Config.h"

#include <cstdlib>
#include <filesystem>
#include <stdexcept>

namespace nuflow {

std::string ResolvePath(const std::string& path) {
  namespace fs = std::filesystem;
  if (fs::exists(path)) return path;
  if (!fs::path(path).is_absolute()) {
    if (const char* dir = std::getenv("NUFLOW_DIR")) {
      fs::path p = fs::path(dir) / path;
      if (fs::exists(p)) return p.string();
    }
    fs::path p = fs::path(NUFLOW_SOURCE_DIR) / path;
    if (fs::exists(p)) return p.string();
  }
  throw std::runtime_error("nuflow: cannot find file '" + path + "' (searched cwd, $NUFLOW_DIR, source dir)");
}

uint64_t HashString(const std::string& s) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (unsigned char c : s) {
    h ^= c;
    h *= 0x100000001b3ULL;
  }
  return h;
}

YAML::Node LoadConfig(const std::string& path) { return YAML::LoadFile(ResolvePath(path)); }

}  // namespace nuflow
