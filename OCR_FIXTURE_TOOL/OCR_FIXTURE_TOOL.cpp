#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
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
bool DEBUG = false;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Mesh::Vertex_index Vertex_index;

bool is_valid_mesh(Mesh mesh) {
	std::stringstream buffer;
	std::streambuf* prevcerr = std::cerr.rdbuf(buffer.rdbuf());
	bool isValid = CGAL::is_valid_polygon_mesh(mesh, DEBUG);
	std::cerr.rdbuf(prevcerr);
	std::string line;
	while (std::getline(buffer, line)) {
		if (DEBUG) std::cout << Yellow << "      Validation output: " << ColorEnd << line << ColorEnd << std::endl;
	}
	return isValid;
}

bool repair_and_validate_mesh(Mesh& mesh) {
	if (DEBUG) std::cout << Yellow << "      Number of removed vertices: " << ColorEnd << PMP::remove_isolated_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Number of new vertices: " << ColorEnd << PMP::duplicate_non_manifold_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Number of stitches: " << ColorEnd << PMP::stitch_borders(mesh) << std::endl;
	return is_valid_mesh(mesh);
}

void get_dimensions(Mesh mesh, double& modelWidth, double& modelLength, double& modelHeight) {
	std::vector<Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
	modelWidth = static_cast<double>(bbox.xmax() - bbox.xmin());
	modelLength = static_cast<double>(bbox.ymax() - bbox.ymin());
	modelHeight = static_cast<double>(bbox.zmax() - bbox.zmin());
	if (DEBUG) std::cout << Yellow << "      Dimensions:" << ColorEnd
		<< "  (W"
		<< modelWidth << "  L" 
		<< modelLength << "  H" 
		<< modelHeight << ")" << std::endl;
}

void get_center(Mesh mesh, Point& center) {
	CGAL::Bbox_3 bbox;
	for (auto v : mesh.vertices()) {
		bbox += mesh.point(v).bbox();
	}
	center = Point((bbox.xmin() + bbox.xmax()) / 2.0, (bbox.ymin() + bbox.ymax()) / 2.0, (bbox.zmin() + bbox.zmax()) / 2.0);
	if (DEBUG) std::cout << Yellow << "      Center:" << ColorEnd
		<< "  ("
		<< center.x() << ", "
		<< center.y() << ", "
		<< center.z() << ")" << std::endl;
}

void get_centroid(Mesh mesh, Point& centroid) {
	std::vector<Point> vertices;
	for (auto v : mesh.vertices()) {
		vertices.push_back(mesh.point(v));
	}
	centroid = CGAL::centroid(vertices.begin(), vertices.end());
	if (DEBUG) std::cout << Yellow << "      Centroid:" << ColorEnd
		<< "  ("
		<< centroid.x() << ", "
		<< centroid.y() << ", "
		<< centroid.z() << ")" << std::endl;
}

void cut_mesh(Mesh& mesh, double model_height, double max_height) {
	double size = 100.0, height = model_height - max_height;  // Adjust based on the expected size of the tm bounding box
	if (height > 0) {
		Mesh clipper, Result_Mesh;
		Vertex_index v0 = clipper.add_vertex(Point(-size, -size, height));
		Vertex_index v1 = clipper.add_vertex(Point(size, -size, height));
		Vertex_index v2 = clipper.add_vertex(Point(size, size, height));
		Vertex_index v3 = clipper.add_vertex(Point(-size, size, height));
		Vertex_index v4 = clipper.add_vertex(Point(-size, -size, -height * 2));
		Vertex_index v5 = clipper.add_vertex(Point(size, -size, -height * 2));
		Vertex_index v6 = clipper.add_vertex(Point(size, size, -height * 2));
		Vertex_index v7 = clipper.add_vertex(Point(-size, size, -height * 2));
		// Bottom face
		clipper.add_face(v0, v1, v2);
		clipper.add_face(v2, v3, v0);
		// Top face
		clipper.add_face(v4, v6, v5);
		clipper.add_face(v6, v4, v7);
		// Four side faces
		clipper.add_face(v0, v4, v1);
		clipper.add_face(v1, v4, v5);
		clipper.add_face(v1, v5, v2);
		clipper.add_face(v2, v5, v6);
		clipper.add_face(v2, v6, v3);
		clipper.add_face(v3, v6, v7);
		clipper.add_face(v3, v7, v0);
		clipper.add_face(v0, v7, v4);

		//Kernel::Plane_3 plane(0, 0, 1, -height); 
		//PMP::clip(mesh, plane, PMP::parameters::clip_volume(true));
		if (DEBUG) std::cout << Yellow << "      Cutting mesh at Z:  " << ColorEnd << height << std::endl;
		if (!PMP::corefine_and_compute_difference(mesh, clipper, Result_Mesh)) {
			std::cerr << Red << "      Cutting mesh failed." << ColorEnd << std::endl;
		}
		mesh.clear();
		mesh = Result_Mesh;
	}
	else {
		if (DEBUG) std::cout << Yellow << "      No Cutting mesh needed:  " << ColorEnd << height << std::endl;
	}
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
	double min_z = std::numeric_limits<double>::infinity();
	for (auto v : mesh.vertices()) {
		double z = mesh.point(v).z();
		if (z < min_z) {
			min_z = z;
		}
	}
	if (DEBUG) std::cout << Yellow << "      Settling mesh at Z:  " << ColorEnd << -min_z << std::endl;
	for (auto v : mesh.vertices()) {
		Point p = mesh.point(v);
		mesh.point(v) = Point(p.x(), p.y(), p.z() - min_z);
	}

	if (DEBUG) std::cout << Yellow << "      Applying translation:  " << ColorEnd << translation_vector << std::endl;
	for (auto v : mesh.vertices()) {
		mesh.point(v) = mesh.point(v) + translation_vector;
	}
}

bool read_STL(const std::string& filename, Mesh& mesh) {
	fs::path filepath(filename);
	mesh.clear();
	if (DEBUG) std::cout << Yellow << "      Reading STL file:  " << ColorEnd << filepath.filename() << std::endl;
	if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
		std::cerr << Red << "Error: Cannot read the STL file:  " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	return true;
}

bool write_STL(const std::string& filename, const Mesh& mesh) {
	fs::path filepath(filename);
	if (DEBUG) std::cout << Yellow << "      Writting STL file:  " << ColorEnd << filepath.filename() << std::endl;
	if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(10))) {
		std::cerr << Red << "Error: Cannot write the STL file:  " << ColorEnd << filepath.filename() << std::endl;
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
		if (std::string(argv[i]) == "-DB") {
			DEBUG = true;
			continue;
		}
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
		std::cerr << Yellow << "Usage: OCR_FIXTURE_TOOL.exe -O out.stl -N id -D Depth [-I model.stl] [-MH MaxHeight] [-DB Debug]   (V3.0 CreatedByBanna)" << ColorEnd << std::endl;
		return EXIT_FAILURE;
	}

	std::string Output_Path_Str = args["-O"], ID_Str = args["-N"], Model_Path_Str = args["-I"], cutting_height_Str = args["-MH"];
	zDepth = 4.0f + std::atof(args["-D"].c_str());

	std::transform(ID_Str.begin(), ID_Str.end(), ID_Str.begin(), [](unsigned char c) { return std::toupper(c); });

	Mesh Letter_Mesh, Tag_Mesh, Fixture_Mesh, Model_Mesh, Result_Mesh;

	for (char c : ID_Str) {
		double FontWidth = 0.0f, FontLength = 0.0f, FontHeight = 0.0f;
		std::string Letter_File = "models/font/" + std::string(1, c) + ".stl";

		if (!read_STL(Letter_File, Letter_Mesh)) return EXIT_FAILURE;

		get_dimensions(Letter_Mesh, FontWidth, FontLength, FontHeight);
		Point centroid, center;
		get_centroid(Letter_Mesh, centroid);
		get_center(Letter_Mesh, center);

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


	if (!Model_Path_Str.empty()) {
		if (!read_STL(Model_Path_Str, Model_Mesh)) return EXIT_FAILURE;
		if (DEBUG) std::cout << Yellow << "      Number of removed vertices: " << ColorEnd << PMP::remove_isolated_vertices(Model_Mesh) << std::endl;

		double Width, Length, Height;
		Point centroid, center;
		get_dimensions(Model_Mesh, Width, Length, Height);
		if (!cutting_height_Str.empty()) {
			cut_mesh(Model_Mesh, Height, std::atof(cutting_height_Str.c_str()));
		}
		//PMP::orient(Model_Mesh);
		get_centroid(Model_Mesh, centroid);
		get_center(Model_Mesh, center);
		translate_mesh(Model_Mesh, Kernel::Vector_3(-centroid.x(), -centroid.y() + 7, 0));

		Mesh Sub_Mesh = Result_Mesh;
		Result_Mesh.clear();

		//Result_Mesh = Model_Mesh;
		if (!PMP::corefine_and_compute_union(Model_Mesh, Sub_Mesh, Result_Mesh)) {
		    std::cerr << Red << "      Model Addition failed." << ColorEnd << std::endl;
			Result_Mesh.clear();
			CGAL::copy_face_graph(Sub_Mesh, Result_Mesh);
			CGAL::copy_face_graph(Model_Mesh, Result_Mesh);
		}
	}

	if (!is_valid_mesh(Result_Mesh)) {
		std::cerr << Red << "      Mesh is not valid." << ColorEnd << std::endl;
		if (repair_and_validate_mesh(Result_Mesh)) {
			if (DEBUG) std::cout << Red << "      Mesh repaired." << ColorEnd << std::endl;
		}
		else {
			std::cerr << Red << "      Failed to repair or validate the mesh." << ColorEnd << std::endl;
			return EXIT_FAILURE;
		}
	}

	if (!write_STL(Output_Path_Str, Result_Mesh)) return EXIT_FAILURE;

	std::cout << Green << "      Operation completed successfully." << ColorEnd << std::endl;
	return EXIT_SUCCESS;
}
