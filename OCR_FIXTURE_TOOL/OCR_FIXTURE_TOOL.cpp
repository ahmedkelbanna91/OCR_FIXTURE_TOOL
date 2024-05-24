#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/bounding_box.h>
#include <CGAL/centroid.h>
#include <algorithm> 
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <filesystem>
#include "include\rang.hpp"

//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };


namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = std::filesystem;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

void get_model_dimensions(Mesh mesh, double& modelWidth, double& modelLength, double& modelHeight) {
	std::vector<Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
	modelWidth = static_cast<double>(bbox.xmax() - bbox.xmin());
	modelLength = static_cast<double>(bbox.ymax() - bbox.ymin());
	modelHeight = static_cast<double>(bbox.zmax() - bbox.zmin());
}

void compute_centroid(Mesh mesh, Point& centroid) {
	std::vector<Point> vertices;
	for (auto v : mesh.vertices()) {
		vertices.push_back(mesh.point(v));
	}
	centroid = CGAL::centroid(vertices.begin(), vertices.end());
}

void scaleMesh(Mesh& mesh, double XYscale, double Zscale, double zThreshold, double XYtopscale) {
	for (auto v : mesh.vertices()) {
		Point& p = mesh.point(v);
		double new_x, new_y, new_z;

		if (p.z() > zThreshold) {
			new_x = p.x() * XYtopscale;
			new_y = p.y() * XYtopscale;
		}
		else {
			new_x = p.x() * XYscale;
			new_y = p.y() * XYscale;
		}
		new_z = p.z() * Zscale;

		mesh.point(v) = Point(new_x, new_y, new_z);
	}
}

void translate_mesh(Mesh& mesh, const Vector& translation_vector) {
	for (auto v : mesh.vertices()) {
		mesh.point(v) = mesh.point(v) + translation_vector;
	}
}

bool read_STL(const std::string& filename, Mesh& mesh) {
	fs::path filepath(filename);
	mesh.clear();
	if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
		std::cerr << Red << "Error: Cannot read the STL file " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	return true;
}

bool write_STL(const std::string& filename, const Mesh& mesh) {
	fs::path filepath(filename);
	if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(10))) {
		std::cerr << Red << "Error: Cannot write the STL file: " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	return true;
}

int main(int argc, char* argv[]) {
	bool lastWasDigit = false;
	double offsetX = 0.0f, offsetY = 0.0f, offsetZ = 0.0f;
	double XYscale = 0.18f, XYtopscale = 0.18f, Zscale = 0.30f;
	double zThreshold = 0.1f;
	double Xspacing = 0.8f, Yspacing = 2.9f;
	double Xtranslate = -6.5f, Ytranslate = -7.5f, zDepth= 4.0f; //2.6f

	std::map<std::string, std::string> args;
	for (int i = 1; i < argc; ++i) {
		if (i + 1 < argc) {
			args[argv[i]] = argv[i + 1];
			i++;
		}
		else {
			std::cerr << Red << "      Missing value for " << ColorEnd << argv[i] << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (args.find("-O") == args.end() || args.find("-N") == args.end() || args.find("-D") == args.end()) {
		std::cerr << Yellow << "Usage: OCR_FIXTURE_TOOL.exe -O out.stl -N id -D Depth [-I model.stl]    (V3.0 CreatedByBanna)" << ColorEnd << std::endl;
		return EXIT_FAILURE;
	}

	std::string Output_Path = args["-O"], ID = args["-N"], Model_Path = args["-I"];
	zDepth = 4.0f + std::atof(args["-D"].c_str());

	std::transform(ID.begin(), ID.end(), ID.begin(), [](unsigned char c) { return std::toupper(c); });

	Mesh Letter_Mesh, Tag_Mesh, Fixture_Mesh, Model_Mesh, Result_Mesh;

	for (char c : ID) {
		double FontWidth = 0.0f, FontLength = 0.0f, FontHeight = 0.0f;
		std::string filename = "models/font/" + std::string(1, c) + ".stl";

		if (!read_STL(filename, Letter_Mesh)) return EXIT_FAILURE;

		get_model_dimensions(Letter_Mesh, FontWidth, FontLength, FontHeight);

		if (std::isdigit(c)) {
			lastWasDigit = true;
		}
		else if (lastWasDigit) {
			offsetY -= (FontLength * XYscale) + Yspacing;
			offsetX = 0.15f;
			lastWasDigit = false;
		}
		scaleMesh(Letter_Mesh, XYscale, Zscale, zThreshold, XYtopscale);
		translate_mesh(Letter_Mesh, Kernel::Vector_3(offsetX, offsetY, offsetZ));
		offsetX += (FontWidth * XYscale) + Xspacing;
		CGAL::copy_face_graph(Letter_Mesh, Tag_Mesh);
	}

	translate_mesh(Tag_Mesh, Kernel::Vector_3(Xtranslate, Ytranslate, zDepth));

	if (!read_STL("models/fixture.stl", Fixture_Mesh)) return EXIT_FAILURE;

	if (!PMP::corefine_and_compute_difference(Fixture_Mesh, Tag_Mesh, Result_Mesh)) {
		std::cerr << Red << "      Subtraction operation failed." << ColorEnd << std::endl;
		return EXIT_FAILURE;
	}


	if (!Model_Path.empty()) {
		if (!read_STL(Model_Path, Model_Mesh)) return EXIT_FAILURE;
		
		Point centroid;
		PMP::orient(Model_Mesh);
		compute_centroid(Model_Mesh, centroid);
		translate_mesh(Model_Mesh, Kernel::Vector_3(-centroid.x(), -centroid.y() + 7, 0));

		Mesh Sub_Mesh = Result_Mesh;
		Result_Mesh.clear();

		if (!PMP::corefine_and_compute_union(Model_Mesh, Sub_Mesh, Result_Mesh)) {
		    std::cerr << Red << "      Model Addition failed." << ColorEnd << std::endl;
			Result_Mesh.clear();
			CGAL::copy_face_graph(Sub_Mesh, Result_Mesh);
			CGAL::copy_face_graph(Model_Mesh, Result_Mesh);
		}
	}

	if (!write_STL(Output_Path, Result_Mesh)) return EXIT_FAILURE;

	std::cout << Green << "      Operation completed successfully." << ColorEnd << std::endl;
	return EXIT_SUCCESS;
}
