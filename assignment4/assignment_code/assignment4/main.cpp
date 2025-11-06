#include <iostream>
#include <chrono>

#include "gloo/Scene.hpp"
#include "gloo/components/MaterialComponent.hpp"

#include "hittable/Sphere.hpp"
#include "Tracer.hpp"
#include "SceneParser.hpp"
#include "ArgParser.hpp"
#include "PerspectiveCamera.hpp"
#include "SpiralCamera.hpp"
#include <memory>

using namespace GLOO;

int main(int argc, const char* argv[]) {
  ArgParser arg_parser(argc, argv);
  SceneParser scene_parser;
  auto scene = scene_parser.ParseScene("assignment4/" + arg_parser.input_file);

  // Create camera based on command-line option
  Camera* camera_ptr;
  std::unique_ptr<Camera> camera_unique_ptr;

  const CameraSpec& camera_spec = scene_parser.GetCameraSpec();

  if (arg_parser.camera_type == "spiral") {
    // Create SpiralCamera
    camera_unique_ptr = std::unique_ptr<Camera>(new SpiralCamera(
        camera_spec.center,
        camera_spec.direction,
        camera_spec.up,
        camera_spec.fov,
        arg_parser.twist));
    camera_ptr = camera_unique_ptr.get();
  } else {
    // Default: create PerspectiveCamera
    camera_unique_ptr = std::unique_ptr<Camera>(new PerspectiveCamera(camera_spec));
    camera_ptr = camera_unique_ptr.get();
  }

  Tracer tracer(camera_ptr,
                glm::ivec2(arg_parser.width, arg_parser.height),
                arg_parser.bounces, scene_parser.GetBackgroundColor(),
                scene_parser.GetCubeMapPtr(), arg_parser.shadows);

  // Set anti-aliasing parameters from command line
  tracer.SetSamplesPerDim(arg_parser.samples_per_dim);

  // Convert string filter type to FilterType enum
  FilterType filter_type = FilterType::Box;  // Default
  if (arg_parser.aa_filter == "tent") {
    filter_type = FilterType::Tent;
  } else if (arg_parser.aa_filter == "box") {
    filter_type = FilterType::Box;
  }
  tracer.SetFilterType(filter_type);

  tracer.Render(*scene, arg_parser.output_file);
  return 0;
}
