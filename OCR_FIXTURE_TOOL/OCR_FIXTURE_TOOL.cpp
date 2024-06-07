#include <future> 
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <limits>
#include <algorithm> 
#include <vector>
#include <string>
#include <cmath>
#include <regex>
#include <map>
#include <windows.h>
#include "rang.hpp"
#include "OCR_font_STL.h"

#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkIdList.h>
#include <vtkAnnotatedCubeActor.h>
#include <vtkSphereSource.h>
#include <vtkAxesActor.h>
#include <vtkDiskSource.h>
#include <vtkCaptionActor2D.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>
#include <vtkHardwarePicker.h>
#include <vtkTextActor.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkSliderWidget.h>
#include <vtkCommand.h>

#include <vtkPlaneWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkButtonWidget.h>
#include <vtkImageData.h>
#include <vtkFreeTypeTools.h>
#include <vtkImageData.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkOpenGLRenderWindow.h>

#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkSplineWidget.h>
#include <vtkSplineWidget2.h>
#include <vtkSplineRepresentation.h>
#include <vtkKochanekSpline.h>



#define M_PI 3.14159265358979323846

//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };

bool DEBUG = true;
namespace fs = std::filesystem;



namespace PMP = CGAL::Polygon_mesh_processing;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Halfedge_index Halfedge_index;
typedef Mesh::Face_index Face_index;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef CGAL::Aff_transformation_3<Kernel> Transformation;

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
class C_InteractorStyle;

vtkNew<vtkActor> staticActor;
vtkNew<vtkActor> movableActor;

double maxcut = 10.0, mincut = 0.0;


class RotationSliderCallback : public vtkCommand {
public:
	static RotationSliderCallback* New() {
		return new RotationSliderCallback;
	}

	RotationSliderCallback() : RotActor(nullptr) {}

	void Execute(vtkObject* caller, unsigned long, void*) override {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		vtkSliderRepresentation* sliderRep = dynamic_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation());
		if (!sliderRep) return;
		double value = scaleValue(sliderRep->GetValue(), 180.0);

		if (RotPtr) *RotPtr = value;  // Update the external CutHeight variable
		else std::cerr << "Warning: RotPtr is not initialized." << std::endl;
		
		if (this->RotActor) {
			double CurrentRot[4];
			this->RotActor->GetOrientation(CurrentRot);
			if (DEBUG) std::cout << Yellow << "      Mesh Orientation:  " << ColorEnd
				<< CurrentRot[0] << "  "
				<< CurrentRot[1] << "  "
				<< CurrentRot[2] << "  "
				<< CurrentRot[3] << std::endl;
			this->RotActor->SetOrientation(CurrentRot[0], CurrentRot[1], value);  // Set Z-position of the cutting plane
			//if (DEBUG) std::cout << Yellow << "      Rotation Slider : " << ColorEnd << value << std::endl;
		}
		char label[50];
		sprintf_s(label, "%.1f", value);  // Format to two decimal places
		sliderRep->SetLabelFormat(label);

		sliderWidget->GetInteractor()->GetRenderWindow()->Render(); // Update the display
	}

	void SetRotActor(vtkActor* actor) {
		this->RotActor = actor;
	}
	void SetRotPtr(double* ptr) {
		RotPtr = ptr;
	}
	double scaleValue(double input, double factor) {
		double normalizedInput = input / factor;  // Normalize to -1 to 1
		return factor * normalizedInput * normalizedInput * (input < 0 ? -1 : 1);  // Scale back to -180 to 180
	}

private:
	vtkActor* RotActor;
	double* RotPtr = nullptr;
};

class CutSliderCallback : public vtkCommand
{
public:
	static CutSliderCallback* New() {
		return new CutSliderCallback();
	}

	CutSliderCallback() : CutHeightPtr(nullptr){}

	virtual void Execute(vtkObject* caller, unsigned long eventId, void*) override {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		vtkSliderRepresentation* sliderRep = dynamic_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation());
		if (!sliderRep) return;
		double value = sliderRep->GetValue();
		if (CutHeightPtr) *CutHeightPtr = value;  // Update the external CutHeight variable
		else std::cerr << "Warning: CutHeightPtr is not initialized." << std::endl;

		if (value >= mincut || value <= maxcut) {
			double CurrentPos[4];
			this->CuttingDisk->GetPosition(CurrentPos);
			this->CuttingDisk->SetPosition(CurrentPos[0], CurrentPos[1], -value);  // Set Z-position of the cutting plane
			//if (DEBUG) std::cout << Yellow << "      Cut Slider : " << ColorEnd << value << std::endl;
		}

		char label[50];
		sprintf_s(label, "%.1f", value);  // Format to two decimal places
		sliderRep->SetLabelFormat(label);

		sliderWidget->GetInteractor()->GetRenderWindow()->Render(); // Update the display
	}

	void SetModelActor(vtkActor* actor) {
		this->CuttingDisk = actor;
	}
	void SetCutHeightPtr(double* ptr) {
		CutHeightPtr = ptr;
	}

private:
	vtkActor* CuttingDisk;
	double* CutHeightPtr = nullptr;
};

class C_InteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
	static C_InteractorStyle* New();
	//vtkTypeMacro(C_InteractorStyle, vtkInteractorStyleTrackballCamera);

	vtkSmartPointer<vtkHardwarePicker> Picker;
	vtkSmartPointer<vtkTextActor> TextActor;
	vtkSmartPointer<vtkActor> SelectedMesh;  // The mesh to pick
	vtkSmartPointer<vtkActor> MeshActor;
	bool IsMeshSelected = false;
	int LastPosition[2] = { -1, -1 };
	double X_offset, Y_offset;

	C_InteractorStyle() : X_offset(0), Y_offset(0) {
		this->Picker = vtkSmartPointer<vtkHardwarePicker>::New();
		this->TextActor = vtkSmartPointer<vtkTextActor>::New();
		this->TextActor->GetTextProperty()->SetFontSize(20);
		this->TextActor->GetTextProperty()->SetColor(0.7, 0.5, 0.3);
		this->TextActor->SetDisplayPosition(10, 10);
		this->TextActor->SetInput("Select Model");

		this->SelectedMesh = nullptr;
		this->IsMeshSelected = false;
		this->LastPosition[0] = this->LastPosition[1] = 0;
	}

	virtual void OnLeftButtonDown() override
	{
		int* clickPos = this->GetInteractor()->GetEventPosition();
		this->Picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
		vtkActor* _actor = vtkActor::SafeDownCast(this->Picker->GetActor());

		if (_actor && _actor == this->MeshActor) {
			this->SelectedMesh = _actor;
			this->IsMeshSelected = true;
			this->LastPosition[0] = clickPos[0];
			this->LastPosition[1] = clickPos[1];

			if (DEBUG) std::cout << Yellow << "      Mesh selected!" << ColorEnd << std::endl;
			this->TextActor->SetInput("Mesh selected!");
			this->SelectedMesh->GetProperty()->SetColor(1.0, 1.0, 0.0);  // Change color when selected
		}
		else {
			if (DEBUG) std::cout << Yellow << "      Nothing selected!" << std::endl;
			this->LastPosition[0] = 0; // Reset last position if no valid selection
			this->LastPosition[1] = 0;
		}
		this->Interactor->GetRenderWindow()->Render();
	}

	virtual void OnLeftButtonUp() override {
		if (this->IsMeshSelected && this->SelectedMesh) {
			double CurrentPos[4];
			this->SelectedMesh->GetPosition(CurrentPos);

			X_offset = CurrentPos[0];
			Y_offset = CurrentPos[1];

			if (DEBUG) std::cout << Yellow << "      Final position offset: " 
				<< ColorEnd << "(" << CurrentPos[0] << ", " << CurrentPos[1] << ")" << std::endl;
			if (DEBUG) std::cout << Yellow << "      Mesh deselected!" << ColorEnd << std::endl;
			this->TextActor->SetInput("Mesh deselected!");
			this->SelectedMesh->GetProperty()->SetColor(0.85, 0.85, 0.85);  // Reset to original color
			this->Interactor->GetRenderWindow()->Render(); // Ensure the scene gets updated

			this->IsMeshSelected = false;
			this->SelectedMesh = nullptr;
		}
	}

	virtual void OnMouseMove() override {
		if (this->IsMeshSelected && this->SelectedMesh) {
			vtkRenderWindowInteractor* rwi = this->Interactor;
			vtkRenderer* renderer = this->GetDefaultRenderer();

			int* newPos = rwi->GetEventPosition();

			// Convert new mouse position to world coordinates
			double displayPos[3] = { double(newPos[0]), double(newPos[1]), 0.0 };
			renderer->SetDisplayPoint(displayPos);
			renderer->DisplayToWorld();
			renderer->GetWorldPoint(displayPos);

			// Convert last position to world coordinates
			double lastDisplayPos[3] = { double(LastPosition[0]), double(LastPosition[1]), 0.0 };
			renderer->SetDisplayPoint(lastDisplayPos);
			renderer->DisplayToWorld();
			renderer->GetWorldPoint(lastDisplayPos);

			// Calculate movement deltas
			double dx = displayPos[0] - lastDisplayPos[0];
			double dy = displayPos[1] - lastDisplayPos[1];

			// Update mesh position
			double pos[3];
			this->SelectedMesh->GetPosition(pos);
			double POSX = pos[0] + dx, POSY = pos[1] + dy;
			this->SelectedMesh->SetPosition(POSX, POSY, pos[2]);

			if (DEBUG) std::cout << Yellow << "      Moved to " << ColorEnd 
				<< "(X " << POSX << ", Y " << POSY << ")" << std::endl;
			this->TextActor->SetInput(("(X " + std::to_string(POSX) + ", Y " + std::to_string(POSY) + ")").c_str());

			this->LastPosition[0] = newPos[0];
			this->LastPosition[1] = newPos[1];

			rwi->Render();
		}

		if (!this->IsMeshSelected) {
			vtkInteractorStyleTrackballCamera::OnMouseMove();
		}
	}

	virtual void OnMiddleButtonDown() override {	
		this->StartPan();
	}

	virtual void OnMiddleButtonUp() override {
		this->EndPan();
	}

	virtual void OnRightButtonDown() override {
		this->StartRotate();
	}

	virtual void OnRightButtonUp() override {
		this->EndRotate();
	}

	virtual void Pan() override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkRenderWindowInteractor* rwi = this->Interactor;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		int* lastPos = rwi->GetLastEventPosition();
		int* newPos = rwi->GetEventPosition();
		double dx = newPos[0] - lastPos[0];
		double dy = newPos[1] - lastPos[1];
		double scale = 0.05; // Adjust this scale to control the sensitivity of panning
		dx *= scale;
		dy *= scale;
		double right[3], up[3];
		camera->GetViewUp(up);
		this->CurrentRenderer->GetActiveCamera()->OrthogonalizeViewUp();
		vtkMath::Cross(camera->GetDirectionOfProjection(), up, right);
		vtkMath::Normalize(right);
		double cameraPosition[3], cameraFocalPoint[3];
		camera->GetPosition(cameraPosition);
		camera->GetFocalPoint(cameraFocalPoint);
		for (int i = 0; i < 3; i++) {
			cameraPosition[i] += dx * right[i] + dy * up[i];
			cameraFocalPoint[i] += dx * right[i] + dy * up[i];
		}
		camera->SetPosition(cameraPosition);
		camera->SetFocalPoint(cameraFocalPoint);
		this->CurrentRenderer->ResetCameraClippingRange();
		rwi->Render();
	}

	virtual void Dolly(double amount) override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		const double DollyScaleFactor = 0.1, minDis = 50, maxDis = 300;
		amount = 1.0 + (amount - 1.0) * DollyScaleFactor;

		double newDistance = (camera->GetDistance()) * (1.0 / amount);  // Adjust the interpretation of amount

		if (newDistance < minDis) camera->SetDistance(minDis); 
		else if (newDistance > maxDis) camera->SetDistance(maxDis);
		else camera->Dolly(amount); 
		
		//if (DEBUG) std::cout << Yellow << "      Current Zoom Level: " << ColorEnd << camera->GetDistance() << std::endl;
		this->CurrentRenderer->ResetCameraClippingRange();
		this->Interactor->Render();
	}
};
vtkStandardNewMacro(C_InteractorStyle);

vtkNew<vtkPolyData> mesh_to_vtk(const Mesh& mesh, bool DEBUG) {
	vtkNew<vtkPoints> points;
	vtkNew<vtkCellArray> polygons;

	for (auto v : mesh.vertices()) {
		const Point& p = mesh.point(v);
		points->InsertNextPoint(p.x(), p.y(), p.z());
	}

	for (auto f : mesh.faces()) {
		vtkNew<vtkIdList> polygon;
		CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(f), mesh);
		for (; vbegin != vend; ++vbegin) {
			polygon->InsertNextId(*vbegin);
		}
		polygons->InsertNextCell(polygon);
	}

	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(points);
	polyData->SetPolys(polygons);
	if (DEBUG) std::cout << Yellow << "      Mesh prepared for Viewer." << ColorEnd << std::endl;
	return polyData;
}


void visualize_mesh(Mesh movableMesh, Mesh staticMesh, double& Xoffset, double& Yoffset,
	double& CutHeight, double& RotZ, bool DEBUG) {

	if (DEBUG) std::cout << Yellow << "      Preparing Mesh Viewer." << ColorEnd << std::endl;
	CutHeight = 0.0; RotZ = 0.0;

	byte WIN_TRANS = 245;
	double bkR = 0.129, bkG = 0.129, bkB = 0.141,
		movR = 0.95, movG = 0.95, movB = 0.95,
		stR = 0.7, stG = 0.5, stB = 0.3;
	//double bkR = 0.2, bkG = 0.2, bkB = 0.3;
	// Main Renderer setup
	vtkNew<vtkRenderer> MainRenderer;
	MainRenderer->SetBackground(bkR, bkG, bkB);
	MainRenderer->SetUseDepthPeeling(1);
	MainRenderer->SetMaximumNumberOfPeels(100);
	MainRenderer->SetOcclusionRatio(0.1);

	// Setup for secondary renderer
	vtkNew<vtkRenderer> InsetRenderer;
	InsetRenderer->SetBackground(bkR, bkG, bkB);
	InsetRenderer->SetUseDepthPeeling(1);
	InsetRenderer->SetMaximumNumberOfPeels(100);
	InsetRenderer->SetOcclusionRatio(0.1);
	InsetRenderer->SetViewport(0.75, 0.0, 1.0, 0.25); // X1 Y1   X2 Y2
	InsetRenderer->EraseOff();

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(700, 600);  
	renderWindow->SetWindowName("ABViewer"); 
	renderWindow->SetAlphaBitPlanes(1);
	renderWindow->SetMultiSamples(0); 
	renderWindow->AddRenderer(MainRenderer);
	renderWindow->AddRenderer(InsetRenderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);


	// Mesh mappers and actors
	vtkNew<vtkPolyDataMapper> staticMapper, movableMapper;
	staticMapper->SetInputData(mesh_to_vtk(staticMesh, DEBUG));
	movableMapper->SetInputData(mesh_to_vtk(movableMesh, DEBUG));

	vtkNew<vtkActor> staticActor, movableActor;
	staticActor->SetMapper(staticMapper);
	movableActor->SetMapper(movableMapper);
	staticActor->GetProperty()->SetColor(stR, stG, stB);
	movableActor->GetProperty()->SetColor(movR, movG, movB);


	MainRenderer->AddActor(staticActor);
	MainRenderer->AddActor(movableActor);
	InsetRenderer->AddActor(movableActor);

	//// Orientation Marker Widget - SetupCubeWidget(renderWindowInteractor, MainRenderer);
	//vtkNew<vtkAnnotatedCubeActor> cubeActor;
	//cubeActor->SetFaceTextScale(0.5);
	//cubeActor->SetXPlusFaceText("X+");
	//cubeActor->SetXMinusFaceText("X-");
	//cubeActor->SetYPlusFaceText("Y+");
	//cubeActor->SetYMinusFaceText("Y-");
	//cubeActor->SetZPlusFaceText("+Z");
	//cubeActor->SetZMinusFaceText("Z-");
	//cubeActor->GetXPlusFaceProperty()->SetColor(1, 0, 0);   // Red for X+
	//cubeActor->GetXMinusFaceProperty()->SetColor(1, 0, 0);  // Red for X-
	//cubeActor->GetYPlusFaceProperty()->SetColor(0, 1, 0);   // Green for Y+
	//cubeActor->GetYMinusFaceProperty()->SetColor(0, 1, 0);  // Green for Y-
	//cubeActor->GetZPlusFaceProperty()->SetColor(0, 0, 1);   // Blue for Z+
	//cubeActor->GetZMinusFaceProperty()->SetColor(0, 0, 1);  // Blue for Z-
	//cubeActor->SetXFaceTextRotation(0);
	//cubeActor->SetYFaceTextRotation(0);
	//cubeActor->SetZFaceTextRotation(-90);
	//cubeActor->GetCubeProperty()->SetColor(0.75, 0.75, 0.75);
	//vtkNew<vtkOrientationMarkerWidget> CubeActorWidget;
	//CubeActorWidget->SetOrientationMarker(cubeActor);
	//CubeActorWidget->SetViewport(0.0, 0.85, 0.15, 1.0);
	//CubeActorWidget->SetInteractor(renderWindowInteractor);
	//CubeActorWidget->SetDefaultRenderer(MainRenderer);
	//CubeActorWidget->SetEnabled(1);
	//CubeActorWidget->InteractiveOff();

	// Center XYZ axes - SetupCenterXYZAxes(renderWindowInteractor, MainRenderer);
	vtkNew<vtkAxesActor> XYZaxes;
	XYZaxes->SetTotalLength(6, 6, 6);
	XYZaxes->SetXAxisLabelText("");
	XYZaxes->SetYAxisLabelText("");
	XYZaxes->SetZAxisLabelText("");
	MainRenderer->AddActor(XYZaxes);

	// Center sphere - SetupSphereActor(MainRenderer);
	vtkNew<vtkSphereSource> sphereSource;
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(0.5);
	sphereSource->Update();
	vtkNew<vtkPolyDataMapper> sphereMapper;
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
	vtkNew<vtkActor> sphereActor;
	sphereActor->SetMapper(sphereMapper);
	sphereActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
	MainRenderer->AddActor(sphereActor);


	// XY Cutting disk - SetupCutdiskActor(InsetRenderer);
	vtkNew<vtkDiskSource> diskSource;
	diskSource->SetInnerRadius(0.0); 
	diskSource->SetOuterRadius(40.0);
	diskSource->SetRadialResolution(50);
	diskSource->SetCircumferentialResolution(50);  
	diskSource->SetNormal(0, 0, 1);
	diskSource->Update();  
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(diskSource->GetOutputPort());
	vtkNew<vtkActor> diskActor;
	diskActor->SetMapper(mapper);
	//diskActor->GetProperty()->SetColor(0.2, 0.5, 0.4);
	diskActor->GetProperty()->SetColor(0.2, 0.2, 0.2);
	diskActor->GetProperty()->SetOpacity(0.5);
	InsetRenderer->AddActor(diskActor);

	// Cutting Slider - SetupCutSliderWidget(renderWindowInteractor, movableActor, CutHeight, mincut, maxcut, mincut);
	vtkNew<vtkSliderRepresentation2D> CutSlider;
	CutSlider->SetMinimumValue(mincut); 
	CutSlider->SetMaximumValue(maxcut); 
	CutSlider->SetValue(mincut);
	CutSlider->GetSliderProperty()->SetColor(0.7, 0.5, 0.3);
	CutSlider->GetTitleProperty()->SetColor(1.0, 1.0, 1.0);
	CutSlider->SetTubeWidth(0.005);
	CutSlider->SetSliderLength(0.05);
	CutSlider->SetSliderWidth(0.02);
	CutSlider->SetEndCapLength(0.005);
	CutSlider->SetEndCapWidth(0.02);
	CutSlider->SetTitleText("Cut");
	CutSlider->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	CutSlider->GetPoint1Coordinate()->SetValue(0.9, 0.20);  // Bottom-left 
	CutSlider->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	CutSlider->GetPoint2Coordinate()->SetValue(0.9, 0.80);  // Top-left
	vtkNew<vtkSliderWidget> CutSliderWidget;
	CutSliderWidget->SetInteractor(renderWindowInteractor);
	CutSliderWidget->SetRepresentation(CutSlider);
	CutSliderWidget->SetAnimationModeToAnimate(); 
	CutSliderWidget->SetEnabled(1);
	vtkNew<CutSliderCallback> CutSlidercallback;
	CutSlidercallback->SetModelActor(movableActor); 
	CutSlidercallback->SetCutHeightPtr(&CutHeight);
	CutSliderWidget->AddObserver(vtkCommand::InteractionEvent, CutSlidercallback);

	// Rotate slider - SetupRotSliderWidget(renderWindowInteractor, movableActor, RotZ, -180, 180, 0);
	vtkNew<vtkSliderRepresentation2D> RotSlider;
	RotSlider->SetMinimumValue(-180); 
	RotSlider->SetMaximumValue(180);
	RotSlider->SetValue(0);
	RotSlider->GetSliderProperty()->SetColor(0.7, 0.5, 0.3); 
	RotSlider->GetTitleProperty()->SetColor(1.0, 1.0, 1.0);
	RotSlider->SetTubeWidth(0.005);
	RotSlider->SetSliderLength(0.05);
	RotSlider->SetSliderWidth(0.02);
	RotSlider->SetEndCapLength(0.005);
	RotSlider->SetEndCapWidth(0.02);
	RotSlider->SetTitleText("Rotate");
	RotSlider->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	RotSlider->GetPoint1Coordinate()->SetValue(0.3, 0.1);
	RotSlider->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	RotSlider->GetPoint2Coordinate()->SetValue(0.7, 0.1);
	vtkNew<vtkSliderWidget> RotSliderWidget;
	RotSliderWidget->SetInteractor(renderWindowInteractor);
	RotSliderWidget->SetRepresentation(RotSlider);
	RotSliderWidget->SetAnimationModeToAnimate();
	RotSliderWidget->SetEnabled(1);
	vtkNew<RotationSliderCallback> RotSlidercallback;
	RotSlidercallback->SetRotActor(movableActor);
	RotSlidercallback->SetRotPtr(&RotZ);
	RotSliderWidget->AddObserver(vtkCommand::InteractionEvent, RotSlidercallback);




	  // Setup the spline widget
	//vtkSmartPointer<vtkSplineWidget2> splineWidget = vtkSmartPointer<vtkSplineWidget2>::New();
	//splineWidget->SetInteractor(renderWindowInteractor);
	//splineWidget->CreateDefaultRepresentation();
	//splineWidget->On();
	


	 // Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(MainRenderer);
	style->MeshActor = movableActor;
	MainRenderer->AddActor(style->TextActor);
	renderWindowInteractor->SetInteractorStyle(style);


	// Camera setups
	MainRenderer->ResetCamera();
	MainRenderer->GetActiveCamera()->SetPosition(0, 0, 160);
	MainRenderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
	MainRenderer->GetActiveCamera()->SetViewUp(0, 1, 0);
	MainRenderer->ResetCameraClippingRange();
	
	InsetRenderer->ResetCamera();
	InsetRenderer->GetActiveCamera()->SetPosition(95, 95, 40);
	InsetRenderer->GetActiveCamera()->SetFocalPoint(0, 0, 20);
	InsetRenderer->GetActiveCamera()->SetViewUp(0, 0, 1);
	InsetRenderer->ResetCameraClippingRange();

	renderWindow->Render();

	// Get the HWND from the VTK Render Window
	HWND hWnd = static_cast<HWND>(renderWindow->GetGenericWindowId());
	SetWindowLong(hWnd, GWL_EXSTYLE, GetWindowLong(hWnd, GWL_EXSTYLE) | WS_EX_LAYERED /* | WS_EX_CLIENTEDGE */);
	SetLayeredWindowAttributes(hWnd, 0, WIN_TRANS, LWA_ALPHA);

	if (DEBUG) std::cout << Yellow << "      Viewer Started." << ColorEnd << std::endl;
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();
	if (DEBUG) std::cout << Yellow << "      Viewer Exited." << ColorEnd << std::endl;

	Xoffset = style->X_offset;
	Yoffset = style->Y_offset;

	if (DEBUG) std::cout << Yellow << "      Offsets: " << ColorEnd << "X" << Xoffset << "  Y" << Yoffset 
		<< "  Z" << CutHeight << "  Rot Z" << RotZ << std::endl;
	std::cout << std::endl;
}


bool is_valid_mesh(Mesh mesh, bool DEBUG) {
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


int close_mesh_hole(Mesh& mesh, bool DEBUG) {
	int holes_closed = 0;
	std::vector<Halfedge_index> border_halfedges;
	for (Halfedge_index h : mesh.halfedges()) {
		if (mesh.is_border(h)) {
			border_halfedges.push_back(h);
		}
	}

	for (Halfedge_index h : border_halfedges) {
		std::vector<Face_index> patch_facets;
		std::vector<Vertex_index> patch_vertices;

		bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
			h,
			CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
			.vertex_output_iterator(std::back_inserter(patch_vertices))));

		if (success) {
			holes_closed++;
		}
		else {
			std::cerr << Red << "      Failed to close mesh holes." << ColorEnd << std::endl;
		}
	}
	return holes_closed;
}

void clean_difference(Mesh& mesh, bool DEBUG) {
	auto vpm = get(CGAL::vertex_point, mesh);
	std::vector<std::size_t> component_ids(num_faces(mesh));
	auto component_map = CGAL::make_property_map(component_ids);
	std::size_t num_components = PMP::connected_components(mesh, component_map, PMP::parameters::vertex_point_map(vpm));

	std::vector<std::size_t> component_sizes(num_components, 0);
	for (face_descriptor fd : faces(mesh)) {
		++component_sizes[component_map[fd]];
	}

	std::size_t largest_component_id = std::distance(component_sizes.begin(), 
		std::max_element(component_sizes.begin(), component_sizes.end()));

	std::vector<Face_index> faces_to_remove;
	for (face_descriptor fd : faces(mesh)) {
		if (component_map[fd] != largest_component_id) {
			faces_to_remove.push_back(fd);
		}
	}

	for (Face_index f : faces_to_remove) {
		mesh.remove_face(f);
	}

	mesh.collect_garbage();
	CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
}

bool repair_and_validate_mesh(Mesh & mesh, bool DEBUG) {
	if (DEBUG) std::cout << Yellow << "      Mesh removed vertices: " 
		<< ColorEnd << PMP::remove_isolated_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Mesh new vertices: " 
		<< ColorEnd << PMP::duplicate_non_manifold_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Mesh stitches: " 
		<< ColorEnd << PMP::stitch_borders(mesh) << std::endl;
	mesh.collect_garbage();
	return is_valid_mesh(mesh, DEBUG);
}

void get_dimensions(Mesh mesh, double& modelWidth, double& modelLength, double& modelHeight, bool DEBUG) {
	std::vector<Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
	modelWidth = static_cast<double>(bbox.xmax() - bbox.xmin());
	modelLength = static_cast<double>(bbox.ymax() - bbox.ymin());
	modelHeight = static_cast<double>(bbox.zmax() - bbox.zmin());
	if (DEBUG) std::cout << Yellow << "      Mesh Dimensions:  W " << ColorEnd
		<< modelWidth << "   L " 
		<< modelLength << "   H " 
		<< modelHeight << std::endl;
}

void get_center(Mesh mesh, Point& center, bool DEBUG) {
	CGAL::Bbox_3 bbox;
	for (auto v : mesh.vertices()) {
		bbox += mesh.point(v).bbox();
	}
	center = Point((bbox.xmin() + bbox.xmax()) / 2.0, (bbox.ymin() + bbox.ymax()) / 2.0, (bbox.zmin() + bbox.zmax()) / 2.0);
	if (DEBUG) std::cout << Yellow << "      Mesh Center:  " << ColorEnd
		<< center.x() << "  "
		<< center.y() << "  "
		<< center.z() << std::endl;
}

void get_centroid(Mesh mesh, Point& centroid, bool DEBUG) {
	std::vector<Point> vertices;
	for (auto v : mesh.vertices()) {
		vertices.push_back(mesh.point(v));
	}
	centroid = CGAL::centroid(vertices.begin(), vertices.end());
	if (DEBUG) std::cout << Yellow << "      Mesh Centroid:  " << ColorEnd
		<< centroid.x() << "  "
		<< centroid.y() << "  "
		<< centroid.z() << std::endl;
}

void settle_mesh_z0(Mesh& mesh, bool DEBUG) {
	double min_z = std::numeric_limits<double>::infinity();
	for (auto v : mesh.vertices()) {
		double z = mesh.point(v).z();
		if (z < min_z) min_z = z;
	}

	Kernel::Vector_3 translation_vector(0, 0, -min_z);
	for (auto v : mesh.vertices()) {
		Point p = mesh.point(v) + translation_vector;
		mesh.point(v) = p;
	}
	if (DEBUG) std::cout << Yellow << "      Mesh Settled at Z:  " << ColorEnd << -min_z << std::endl;
}

void cut_mesh(Mesh& mesh, double model_height, double max_height, bool DEBUG) {
	double size = 100.0, height = model_height - max_height , bottom_z = -10;
	if (height >= 0) {
		Mesh clipper, Result_Mesh;
		Vertex_index v0 = clipper.add_vertex(Point(-size, -size, height));
		Vertex_index v1 = clipper.add_vertex(Point(size, -size, height));
		Vertex_index v2 = clipper.add_vertex(Point(size, size, height));
		Vertex_index v3 = clipper.add_vertex(Point(-size, size, height));
		Vertex_index v4 = clipper.add_vertex(Point(-size, -size, bottom_z));
		Vertex_index v5 = clipper.add_vertex(Point(size, -size, bottom_z));
		Vertex_index v6 = clipper.add_vertex(Point(size, size, bottom_z));
		Vertex_index v7 = clipper.add_vertex(Point(-size, size, bottom_z));
		// Top face
		clipper.add_face(v0, v1, v2);
		clipper.add_face(v2, v3, v0);
		// Bottom face
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

		if (!PMP::corefine_and_compute_difference(mesh, clipper, Result_Mesh)) {
			std::cerr << Red << "      Cutting mesh failed." << ColorEnd << std::endl;
		}
		settle_mesh_z0(Result_Mesh, DEBUG);
		mesh.clear();
		mesh = Result_Mesh;
		if (DEBUG) std::cout << Yellow << "      Mesh Cut at Z:  " << ColorEnd << height << std::endl;
	}
	else {
		if (DEBUG) std::cout << Yellow << "      No Cutting mesh needed:  " << ColorEnd << height << std::endl;
	}
}

void extrude_bottom_faces(Mesh& mesh, double target_z, bool DEBUG) {
	double z_threshold = 0.1;
	for (Vertex_index v : mesh.vertices()) {
		Point& p = mesh.point(v);
		if (p.z() >= z_threshold) {
			mesh.point(v) = Point(p.x(), p.y(), p.z() - target_z);
		}
	}
	if (DEBUG) std::cout << Yellow << "      Mesh Extruded:  " << ColorEnd << target_z << std::endl;
}

void scaleMesh(Mesh& mesh, double XYscale, double XYtopscale, double Zscale, double zThreshold, bool DEBUG) {
	for (auto v : mesh.vertices()) {
		Point& point = mesh.point(v);
		double new_x, new_y, new_z;

		if (point.z() > zThreshold) {
			new_x = point.x() * XYtopscale;
			new_y = point.y() * XYtopscale;
		} else {
			new_x = point.x() * XYscale;
			new_y = point.y() * XYscale;
		}
		new_z = point.z() * Zscale;
		mesh.point(v) = Point(new_x, new_y, new_z);
	}
}

void translate_mesh(Mesh& mesh, const Vector& translation_vector, bool DEBUG) {
	for (auto v : mesh.vertices()) {
		mesh.point(v) = mesh.point(v) + translation_vector;
	}
	if (DEBUG) std::cout << Yellow << "      translation Applied:  " << ColorEnd << translation_vector << std::endl;
}

void rotate_mesh(Mesh& mesh, double x_deg, double y_deg, double z_deg, bool DEBUG) {
	double rot_x = x_deg * M_PI / 180.0;
	double rot_y = y_deg * M_PI / 180.0;
	double rot_z = z_deg * M_PI / 180.0;

	double cos_x = std::cos(rot_x), sin_x = std::sin(rot_x);
	Transformation rot_mtx_x(1, 0, 0, 0, 0, cos_x, -sin_x, 0, 0, sin_x, cos_x, 0, 1);
	double cos_y = std::cos(rot_y), sin_y = std::sin(rot_y);
	Transformation rot_mtx_y(cos_y, 0, sin_y, 0, 0, 1, 0, 0, -sin_y, 0, cos_y, 0, 1);
	double cos_z = std::cos(rot_z), sin_z = std::sin(rot_z);
	Transformation rot_mtx_z(cos_z, -sin_z, 0, 0, sin_z, cos_z, 0, 0, 0, 0, 1, 0, 1);

	Transformation combined = rot_mtx_x * rot_mtx_y * rot_mtx_z;
	PMP::transform(combined, mesh);
	if (DEBUG) std::cout << Yellow << "      Rotation Applied:  "
		<< ColorEnd << "X " << x_deg << ", Y " << y_deg << ", Z " << z_deg << std::endl;
}

bool read_STL_data(const std::string& identifier, Mesh& mesh, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	mesh.clear();
	for (const auto& data : FONT_STL) {
		if (data.key == identifier) {
			std::istringstream iss(std::string(reinterpret_cast<const char*>(data.data), data.size), std::ios::binary);
			if (CGAL::IO::read_STL(iss, mesh)) {
				auto finish_ = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed_ = finish_ - start_;
				if (DEBUG) std::cout << "   >> " << Yellow << "Read Mesh: " << ColorEnd << identifier
					<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
				return true;
			}
			break;
		}
	}
	std::cerr << Red << "      Error: No STL data available for:  " << ColorEnd << identifier << std::endl;
	return false;
}

bool read_STL(const std::string& filename, Mesh& mesh , bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	mesh.clear();
	fs::path filepath(filename);
	if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
		std::cerr << Red << "      Error: Cannot read the STL file:  " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	if (DEBUG) std::cout << "   >> " << Yellow << "Read STL File: " << ColorEnd << filepath.filename() 
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	return true;
}

bool write_STL(const std::string& filename, const Mesh& mesh, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	fs::path filepath(filename);
	if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(10))) {
		std::cerr << Red << "      Error: Cannot write the STL file:  " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	if (DEBUG) std::cout << "   << " << Yellow << "Written STL File: " << ColorEnd << filepath.filename()
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	return true;
}

void create_fixture(std::string ID_Str, Mesh Fix_Mesh, Mesh& Res_Mesh, bool DEBUG) {
	bool lastWasDigit = false;
	double offsetX = -6.5, offsetY = -7.5, offsetZ = 4.0;
	double XYscale = 0.18, XYtopscale = 0.18, Zscale = 0.30;
	double zThreshold = 0.1;
	double Xspacing = 0.8, Yspacing = 2.9;
	double zDepth = -0.7;
	Mesh Tag_Mesh;

	if (DEBUG) std::cout << Yellow << "      Creating tag on Fixture:  " << ColorEnd << ID_Str << std::endl;

	std::transform(ID_Str.begin(), ID_Str.end(), ID_Str.begin(), [](unsigned char c) { return std::toupper(c); });

	for (char c : ID_Str) {
		Mesh Letter_Mesh;
		double FontWidth = 0.0, FontLength = 0.0, FontHeight = 0.0;

		if (!read_STL_data(std::string(1, c), Letter_Mesh, false)) continue;

		get_dimensions(Letter_Mesh, FontWidth, FontLength, FontHeight, false);

		if (std::isdigit(c)) {
			lastWasDigit = true;
		} else if (lastWasDigit) {
			offsetY -= (FontLength * XYscale) + Yspacing;
			offsetX = -6.35; // 0.15
			lastWasDigit = false;
		}

		scaleMesh(Letter_Mesh, XYscale, XYtopscale, Zscale, zThreshold, false);
		translate_mesh(Letter_Mesh, Kernel::Vector_3(offsetX, offsetY, offsetZ + zDepth), false);
		offsetX += (FontWidth * XYscale) + Xspacing;
		CGAL::copy_face_graph(Letter_Mesh, Tag_Mesh);
	}

	Res_Mesh.clear();
	if (!PMP::corefine_and_compute_difference(Fix_Mesh, Tag_Mesh, Res_Mesh)) {
		std::cerr << Red << "      Tag Subtraction failed for: " << ColorEnd << ID_Str << std::endl;
	}
}



bool modify_model_mesh(Mesh M_Mesh, Mesh F_Mesh, std::string O_Path, std::string ID_, 
	double M_Xoffset, double M_Yoffset, double C_height, double M_Zrot, bool DEBUG) {

	auto start_ = std::chrono::high_resolution_clock::now();
	Mesh F_ID_Mesh, Result_Mesh;

	std::cout << Yellow << "      Modifying " << ColorEnd << ID_ << std::endl;

	create_fixture(ID_, F_Mesh, F_ID_Mesh, DEBUG);

	if (M_Zrot < -0.001 || M_Zrot > 0.001)
		rotate_mesh(M_Mesh, 0, 0, M_Zrot, DEBUG);

	Point centroid;
	get_centroid(M_Mesh, centroid, DEBUG);
	translate_mesh(M_Mesh, Kernel::Vector_3(-centroid.x() + M_Xoffset, -centroid.y() + M_Yoffset, 0), DEBUG);

	if (C_height < -0.0001 || C_height > 0.0001)
		cut_mesh(M_Mesh, C_height, 0, DEBUG);
	else
		settle_mesh_z0(M_Mesh, DEBUG);
	

	if (!PMP::corefine_and_compute_union(M_Mesh, F_ID_Mesh, Result_Mesh)) {
		std::cerr << Red << "      Model Addition to fixture failed for: " << ColorEnd << ID_ << std::endl;
		Result_Mesh.clear();
		CGAL::copy_face_graph(F_ID_Mesh, Result_Mesh);
		CGAL::copy_face_graph(M_Mesh, Result_Mesh);
	}


	//if (repair_and_validate_mesh(Result_Mesh, DEBUG)) 
	//	if (DEBUG) std::cout << Green << "      Mesh repaired." << ColorEnd << std::endl;
	//else std::cerr << Red << "      Failed to repair or validate the mesh." << ColorEnd << std::endl;


	if (!write_STL(O_Path, Result_Mesh, DEBUG)) return false;

	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	std::cout << " >><< " << Green << "Operation completed successfully." << ColorEnd
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	std::cout << std::endl;
	return true;
}


void promptForNumbers(const std::string& prompt, int& outValue) {
	std::string line;
	bool valid = false;

	while (!valid) {
		std::cout << Yellow << prompt << ColorEnd;
		std::getline(std::cin, line);
		if (line.empty()) {
			outValue = 0;
			valid = true;
		} else {
			size_t dbIndex = line.find("-DB");
			if (dbIndex != std::string::npos) {
				std::string numericPart = line.substr(0, dbIndex);
				std::stringstream numCheck(numericPart);
				int num;
				if (numCheck >> num && numCheck.eof()) {
					DEBUG = true;
					std::cout << Cyan << "      Debug Enabled" << ColorEnd << std::endl;
					outValue = num;
					valid = true;
					continue;
				}
			}
			
			std::stringstream ss(line);
			if (ss >> outValue && ss.eof()) {
				valid = true;
			}
			else {
				std::cout << Red << "      Invalid input:  " << ColorEnd << line << std::endl;
				ss.clear();
			}
		}
	}
}

static void displayUserName() {
	char* username = nullptr;
	char* userdomain = nullptr;
	size_t sizeUsername = 0;
	size_t sizeUserdomain = 0;

	_dupenv_s(&username, &sizeUsername, "USERNAME");
	_dupenv_s(&userdomain, &sizeUserdomain, "USERDOMAIN");

	std::cout << "\n      USERNAME: "
		<< (userdomain ? userdomain : "Unknown") << "\\"
		<< (username ? username : "Unknown") << std::endl;

	free(username);
	free(userdomain);
}


void setConsoleSize(int width, int height) {
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); // Get the standard output handle

	COORD newSize;
	newSize.X = width;
	newSize.Y = 32766; // Maximal possible height for the console window
	SetConsoleScreenBufferSize(hStdout, newSize);

	SMALL_RECT windowSize;
	windowSize.Top = 0;
	windowSize.Left = 0;
	windowSize.Right = width - 1;  // Width of the window
	windowSize.Bottom = height - 1;  // Height of the window

	if (!SetConsoleWindowInfo(hStdout, TRUE, &windowSize)) {
		std::cerr << "Setting console window size failed." << std::endl;
	}
}

int main(int argc, char* argv[]) {
	setConsoleSize(73, 35);

	std::cout << Cyan << "\n===========================" << ColorEnd
		<< Yellow <<	 "'Created by Banna'" << ColorEnd
		<< Cyan <<		 "===========================" << ColorEnd << std::endl;

	std::cout << Cyan << "======================" << ColorEnd
		<< Yellow <<	 "'AB FIXTURE CREATOR TOOL V3'" << ColorEnd
		<< Cyan <<		 "======================" << ColorEnd << std::endl;
	std::cout << Cyan << "========================================================================\n" << ColorEnd << std::endl;


	std::string Main_path = fs::current_path().string();
	std::string outputPath = Main_path + "/output/";
	std::string inputPath = Main_path + "/input/";

	if (!fs::exists(inputPath)) fs::create_directory(inputPath);

	if (!fs::exists(outputPath)) {
		fs::create_directory(outputPath);
	}
	else {
		for (const auto& entry : fs::directory_iterator(outputPath))
			fs::remove_all(entry.path());
	}


	//int caseID;
	//promptForNumbers("      What is the Case ID? ", caseID);
	//
	//std::regex filenameRegex("([0-9]+)(LN|LP|LR|LT|UN|UP|UR|UT)([0-9]+)"); // Regex to split digits, letters, digits
	//for (const auto& entry : fs::directory_iterator(inputPath)) {
	//	if (entry.path().extension() == ".stl") {
	//		std::string filename_WoE = entry.path().stem().string(); // Get filename without extension
	//		std::string filename_WE = entry.path().filename().string(); // Get filename with extension
	//		
	//		std::smatch matches;
	//		if (std::regex_match(filename_WoE, matches, filenameRegex)) {
	//			if (matches.size() == 4) { // Full match + four groups
	//				std::string newFilename = std::to_string(caseID) + matches[2].str() + matches[3].str() + ".stl";
	//
	//				fs::rename(entry.path(), entry.path().parent_path() / newFilename);
	//				std::cout << "   >> " << Yellow << filename_WE << ColorEnd
	//					<< " Renamed to " << Green << newFilename << ColorEnd << std::endl;
	//			}
	//		}
	//		else {
	//			std::cout << Red << "      Nothing matches" << ColorEnd << "  'LN, LP, LR, LT, UN, UP, UR, UT'" << std::endl;
	//		}
	//	}
	//}
	//std::cout << std::endl;
	
	std::vector<std::string> filenames;
	for (const auto& entry : fs::directory_iterator(inputPath)) {
		if (entry.path().extension() == ".stl") { 
			std::string filename_WE= entry.path().stem().string(); 
			filenames.push_back(filename_WE);
		}
	}
	

	//
    //
	// modify here to open the first model in the list in general, also make the rotation and xy translation with vtk realtime render
    //
	//

	auto start = std::chrono::high_resolution_clock::now();

	Mesh Fixture_Mesh;
	double Model_Xoffset = 0.0, Model_Yoffset = 0.0, cut_height = 0.0, Model_Zrot = 0.0;
	if (!read_STL_data("fixture", Fixture_Mesh, DEBUG)) return EXIT_FAILURE;

	// LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER
	std::cout << std::endl;
	std::vector<std::string> filtered_Lower;
	std::copy_if(filenames.begin(), filenames.end(), std::back_inserter(filtered_Lower),
		[](const std::string& name) {
			return name.find("L") != std::string::npos;
		});
	std::sort(filtered_Lower.begin(), filtered_Lower.end());
	for (const auto& LOWER: filtered_Lower) {
		Mesh Model_Mesh, Visual_mesh;
		std::string Model_In_Path = inputPath + LOWER + ".stl";
		std::string Model_Out_Path = outputPath + LOWER + ".stl"; 

		if (!read_STL(Model_In_Path, Model_Mesh, DEBUG)) return EXIT_FAILURE;
		Visual_mesh = Model_Mesh;

		if (std::stoi(LOWER.substr(LOWER.find("LN") + 2)) == 1) { 
			Point centroid;
			get_centroid(Visual_mesh, centroid, DEBUG);
			translate_mesh(Visual_mesh, Kernel::Vector_3(-centroid.x(), -centroid.y(), 0), DEBUG);
			visualize_mesh(Visual_mesh, Fixture_Mesh, Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, true);
		}

		modify_model_mesh(Model_Mesh, Fixture_Mesh, Model_Out_Path, LOWER, 
			Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, DEBUG);
	}

	// UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER 
	std::cout << std::endl;
	std::vector<std::string> filtered_Upper;
	std::copy_if(filenames.begin(), filenames.end(), std::back_inserter(filtered_Upper), 
		[](const std::string& name) {
			return name.find("U") != std::string::npos;
		});
	std::sort(filtered_Upper.begin(), filtered_Upper.end());
	for (const auto& UPPER : filtered_Upper) {
		Mesh Model_Mesh, Visual_mesh;
		std::string Model_In_Path = inputPath + UPPER + ".stl";
		std::string Model_Out_Path = outputPath + UPPER + ".stl";
		
		if (!read_STL(Model_In_Path, Model_Mesh, DEBUG)) return EXIT_FAILURE;
		Visual_mesh = Model_Mesh;

		if (std::stoi(UPPER.substr(UPPER.find("UN") + 2)) == 1) {
			Point centroid;
			get_centroid(Visual_mesh, centroid, DEBUG);
			translate_mesh(Visual_mesh, Kernel::Vector_3(-centroid.x(), -centroid.y(), 0), DEBUG);
			visualize_mesh(Visual_mesh, Fixture_Mesh, Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, true);
		}

		modify_model_mesh(Model_Mesh, Fixture_Mesh, Model_Out_Path, UPPER, 
			Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, DEBUG);
	}
	
	std::cout << std::endl;

	displayUserName();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << Yellow << "      Elapsed time: " << elapsed.count() << " seconds" << ColorEnd << std::endl;

	std::cout << std::endl;
	std::cout << "      Press enter to continue...";
	std::cin.get();  // Waits for the user to press Enter
	return EXIT_SUCCESS;
}

