#pragma once
// YAML configuration helpers and data-path resolution.

#include <cstdint>
#include <string>

#include <yaml-cpp/yaml.h>

namespace nuflow {

// Resolve a (possibly relative) data path: tries the path as given, then $NUFLOW_DIR/path,
// then the source directory the library was built from. Throws if none exists.
std::string ResolvePath(const std::string& path);

// Stable 64-bit FNV-1a hash of a string (used to key cross-section caches).
uint64_t HashString(const std::string& s);

YAML::Node LoadConfig(const std::string& path);

template <typename T>
T Get(const YAML::Node& node, const std::string& key, const T& fallback) {
  if (node && node[key]) return node[key].as<T>();
  return fallback;
}

}  // namespace nuflow
