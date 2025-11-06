#ifndef ARG_PARSER_H_
#define ARG_PARSER_H_

#include <string>

class ArgParser {
 public:
  ArgParser(int argc, const char* argv[]);

  std::string input_file;
  std::string output_file;
  std::string depth_file;
  std::string normals_file;
  size_t width;
  size_t height;

  // Rendering options.
  float depth_min;
  float depth_max;
  size_t bounces;
  bool shadows;

  // Supersampling/anti-aliasing.
  bool jitter;
  bool filter;
  int samples_per_dim;  // Samples per dimension (samples_per_dim^2 total samples per pixel)
  std::string aa_filter;  // Anti-aliasing filter type: "box" or "tent"

  // Camera options.
  std::string camera_type;  // Camera type: "perspective" (default) or "spiral"
  float twist;  // Twist strength for spiral camera (in radians)

 private:
  void SetDefaultValues();
};

#endif  // ARG_PARSER_H
