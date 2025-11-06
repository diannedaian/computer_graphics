#include "ArgParser.hpp"

#include <cstring>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>

ArgParser::ArgParser(int argc, const char* argv[]) {
  SetDefaultValues();

  for (int i = 1; i < argc; i++) {
    // rendering output
    if (!strcmp(argv[i], "-input")) {
      i++;
      assert(i < argc);
      input_file = argv[i];
    } else if (!strcmp(argv[i], "-output")) {
      i++;
      assert(i < argc);
      output_file = argv[i];
    } else if (!strcmp(argv[i], "-size")) {
      i++;
      assert(i < argc);
      width = atoi(argv[i]);
      i++;
      assert(i < argc);
      height = atoi(argv[i]);
    } else if (!strcmp(argv[i], "-bounces")) {
      i++;
      assert(i < argc);
      bounces = atoi(argv[i]);
    } else if (!strcmp(argv[i], "-shadows")) {
      shadows = true;
    } else if (!strcmp(argv[i], "-aa") || !strcmp(argv[i], "-spp")) {
      i++;
      assert(i < argc);
      samples_per_dim = atoi(argv[i]);
      // If N <= 1, treat as no AA (1 sample per pixel)
      if (samples_per_dim <= 1) {
        samples_per_dim = 1;
      }
    } else if (!strcmp(argv[i], "-aafilter")) {
      i++;
      assert(i < argc);
      aa_filter = argv[i];
      // Validate filter type
      if (aa_filter != "box" && aa_filter != "tent") {
        printf("Unknown filter type: '%s'. Use 'box' or 'tent'.\n", aa_filter.c_str());
        exit(1);
      }
    } else if (!strcmp(argv[i], "-camera")) {
      i++;
      assert(i < argc);
      camera_type = argv[i];
      // Validate camera type
      if (camera_type != "perspective" && camera_type != "spiral") {
        printf("Unknown camera type: '%s'. Use 'perspective' or 'spiral'.\n", camera_type.c_str());
        exit(1);
      }
    } else if (!strcmp(argv[i], "-twist")) {
      i++;
      assert(i < argc);
      twist = static_cast<float>(atof(argv[i]));
    } else {
      printf("Unknown command line argument %d: '%s'\n", i, argv[i]);
      exit(1);
    }
  }

  std::cout << "Args:\n";
  std::cout << "- input: " << input_file << std::endl;
  std::cout << "- output: " << output_file << std::endl;
  std::cout << "- width: " << width << std::endl;
  std::cout << "- height: " << height << std::endl;
  std::cout << "- bounces: " << bounces << std::endl;
  std::cout << "- shadows: " << shadows << std::endl;
  std::cout << "- samples_per_dim: " << samples_per_dim << std::endl;
  std::cout << "- aa_filter: " << aa_filter << std::endl;
  std::cout << "- camera_type: " << camera_type << std::endl;
  std::cout << "- twist: " << twist << std::endl;
}

void ArgParser::SetDefaultValues() {
  input_file = "";
  output_file = "";
  normals_file = "";
  width = 200;
  height = 200;

  bounces = 0;
  shadows = false;

  // Anti-aliasing defaults (no AA by default - matches original behavior)
  samples_per_dim = 1;
  aa_filter = "box";

  // Camera defaults
  camera_type = "perspective";  // Default to perspective camera
  twist = 0.5f;  // Default twist strength for spiral camera (in radians)
}
