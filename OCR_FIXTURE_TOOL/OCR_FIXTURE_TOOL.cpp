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
#include <map>
#include "include\rang.hpp"

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


namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = std::filesystem;
bool DEBUG = false;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
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

vtkNew<vtkPolyDataMapper> staticMapper;
vtkNew<vtkPolyDataMapper> movableMapper;
vtkNew<vtkActor> staticActor;
vtkNew<vtkActor> movableActor;

double maxcut = 8.0, mincut = 0.0;

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
		double value = scaleValue(sliderRep->GetValue());

		if (RotPtr) *RotPtr = value;  // Update the external CutHeight variable
		else std::cerr << "Warning: RotPtr is not initialized." << std::endl;
		
		if (this->RotActor) {
			double CurrentRot[4];
			this->RotActor->GetOrientation(CurrentRot);
			this->RotActor->SetOrientation(CurrentRot[0], CurrentRot[1], value);  // Set Z-position of the cutting plane
			if (DEBUG) std::cout << Yellow << "      Rotation Slider : " << ColorEnd << value << std::endl;
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
	double scaleValue(double input) {
		double normalizedInput = input / 180.0;  // Normalize to -1 to 1
		return 180.0 * normalizedInput * normalizedInput * (input < 0 ? -1 : 1);  // Scale back to -180 to 180
		//return 180.0 * std::pow(normalizedInput, 4) * (input < 0 ? -1 : 1);  // Scale back to -180 to 180 with sign
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
			if (DEBUG) std::cout << Yellow << "      Cut Slider : " << ColorEnd << value << std::endl;
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

			if (DEBUG) std::cout << Yellow << "      Final position offset: " << ColorEnd << "(" << CurrentPos[0] << ", " << CurrentPos[1] << ")" << std::endl;
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

			if (DEBUG) std::cout << Yellow << "      Moved to " << ColorEnd << "(X " << POSX << ", Y " << POSY << ")" << std::endl;
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
		
		if (DEBUG) std::cout << Yellow << "      Current Zoom Level: " << ColorEnd << camera->GetDistance() << std::endl;
		this->CurrentRenderer->ResetCameraClippingRange();
		this->Interactor->Render();
	}
};
vtkStandardNewMacro(C_InteractorStyle);

vtkNew<vtkPolyData> mesh_to_vtk(const Mesh& mesh) {
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


void visualize_mesh(Mesh staticMesh, Mesh movableMesh, double& Xoffset, double& Yoffset, double& CutHeight, double& RotZ) {
	if (DEBUG) std::cout << Yellow << "      Preparing Mesh Viewer." << ColorEnd << std::endl;
	CutHeight = 0.0; RotZ = 0.0;

	// Main Renderer setup
	vtkNew<vtkRenderer> MainRenderer;
	MainRenderer->SetBackground(0.2, 0.2, 0.3);
	MainRenderer->SetUseDepthPeeling(1);
	MainRenderer->SetMaximumNumberOfPeels(100);
	MainRenderer->SetOcclusionRatio(0.1);

	// Setup for secondary renderer
	vtkNew<vtkRenderer> insetRenderer;
	insetRenderer->SetBackground(0.2, 0.2, 0.3); 
	insetRenderer->SetUseDepthPeeling(1);
	insetRenderer->SetMaximumNumberOfPeels(100);
	insetRenderer->SetOcclusionRatio(0.1);
	insetRenderer->SetViewport(0.75, 0.0, 1.0, 0.25); // X1 Y1   X2 Y2
	insetRenderer->EraseOff();

	// Camera setups
	//vtkNew<vtkCamera> camera;
	//camera->SetPosition(0, 0, 150);
	//camera->SetFocalPoint(0, 0, 0);
	/*camera->SetViewUp(0, 1, 0);
	MainRenderer->SetActiveCamera(camera);*/

	vtkNew<vtkCamera> sideCamera;
	sideCamera->SetPosition(95, 95, 40);
	sideCamera->SetFocalPoint(0, 0, 20);
	sideCamera->SetViewUp(0, 0, 1);
	insetRenderer->SetActiveCamera(sideCamera);

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(700, 600);  
	renderWindow->SetWindowName("ABViewer"); 
	renderWindow->SetAlphaBitPlanes(1);
	renderWindow->SetMultiSamples(0); 
	renderWindow->AddRenderer(MainRenderer);
	renderWindow->AddRenderer(insetRenderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);


	// Mesh mappers and actors
	vtkNew<vtkPolyDataMapper> staticMapper, movableMapper;
	staticMapper->SetInputData(mesh_to_vtk(staticMesh));
	movableMapper->SetInputData(mesh_to_vtk(movableMesh));

	vtkNew<vtkActor> staticActor, movableActor;
	staticActor->SetMapper(staticMapper);
	movableActor->SetMapper(movableMapper);
	staticActor->GetProperty()->SetColor(0.7, 0.5, 0.3);
	movableActor->GetProperty()->SetColor(0.85, 0.85, 0.85);

	MainRenderer->AddActor(staticActor);
	MainRenderer->AddActor(movableActor);
	insetRenderer->AddActor(staticActor);
	insetRenderer->AddActor(movableActor);


	// Orientation Marker Widget - SetupCubeWidget(renderWindowInteractor, MainRenderer);
	vtkNew<vtkAnnotatedCubeActor> cubeActor;
	cubeActor->SetFaceTextScale(0.5);
	cubeActor->SetXPlusFaceText("X+");
	cubeActor->SetXMinusFaceText("X-");
	cubeActor->SetYPlusFaceText("Y+");
	cubeActor->SetYMinusFaceText("Y-");
	cubeActor->SetZPlusFaceText("+Z");
	cubeActor->SetZMinusFaceText("Z-");
	cubeActor->GetXPlusFaceProperty()->SetColor(1, 0, 0);   // Red for X+
	cubeActor->GetXMinusFaceProperty()->SetColor(1, 0, 0);  // Red for X-
	cubeActor->GetYPlusFaceProperty()->SetColor(0, 1, 0);   // Green for Y+
	cubeActor->GetYMinusFaceProperty()->SetColor(0, 1, 0);  // Green for Y-
	cubeActor->GetZPlusFaceProperty()->SetColor(0, 0, 1);   // Blue for Z+
	cubeActor->GetZMinusFaceProperty()->SetColor(0, 0, 1);  // Blue for Z-
	cubeActor->SetXFaceTextRotation(0);
	cubeActor->SetYFaceTextRotation(0);
	cubeActor->SetZFaceTextRotation(-90);
	cubeActor->GetCubeProperty()->SetColor(0.75, 0.75, 0.75);
	vtkNew<vtkOrientationMarkerWidget> CubeActorWidget;
	CubeActorWidget->SetOrientationMarker(cubeActor);
	CubeActorWidget->SetViewport(0.0, 0.85, 0.15, 1.0);
	CubeActorWidget->SetInteractor(renderWindowInteractor);
	CubeActorWidget->SetDefaultRenderer(MainRenderer);
	CubeActorWidget->SetEnabled(1);
	CubeActorWidget->InteractiveOff();

	// Center XYZ axes - SetupCenterXYZAxes(renderWindowInteractor, MainRenderer);
	vtkNew<vtkAxesActor> XYZaxes;
	XYZaxes->SetTotalLength(10, 10, 10);
	XYZaxes->GetXAxisCaptionActor2D()->SetVisibility(0);
	XYZaxes->GetYAxisCaptionActor2D()->SetVisibility(0);
	XYZaxes->GetZAxisCaptionActor2D()->SetVisibility(0);
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


	// XY Cutting disk - SetupCutdiskActor(insetRenderer);
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
	diskActor->GetProperty()->SetColor(0.2, 0.5, 0.4);
	diskActor->GetProperty()->SetOpacity(0.5);
	insetRenderer->AddActor(diskActor);

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
	CutSlider->GetPoint1Coordinate()->SetValue(0.07, 0.20);  // Bottom-left 
	CutSlider->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	CutSlider->GetPoint2Coordinate()->SetValue(0.07, 0.80);  // Top-left
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


	 // Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(MainRenderer);
	style->MeshActor = movableActor;
	MainRenderer->AddActor(style->TextActor);
	renderWindowInteractor->SetInteractorStyle(style);
	

	MainRenderer->ResetCamera();
	renderWindow->Render();
	renderWindowInteractor->Initialize();

	if (DEBUG) std::cout << Yellow << "      Viewer Started." << ColorEnd << std::endl;
	renderWindowInteractor->Start();
	if (DEBUG) std::cout << Yellow << "      Viewer Exited." << ColorEnd << std::endl;

	Xoffset = style->X_offset;
	Yoffset = style->Y_offset;
}

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


int close_mesh_hole(Mesh& mesh) {
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

void clean_difference(Mesh& mesh) {
	auto vpm = get(CGAL::vertex_point, mesh);
	std::vector<std::size_t> component_ids(num_faces(mesh));
	auto component_map = CGAL::make_property_map(component_ids);
	std::size_t num_components = PMP::connected_components(mesh, component_map, PMP::parameters::vertex_point_map(vpm));

	std::vector<std::size_t> component_sizes(num_components, 0);
	for (face_descriptor fd : faces(mesh)) {
		++component_sizes[component_map[fd]];
	}

	std::size_t largest_component_id = std::distance(component_sizes.begin(), std::max_element(component_sizes.begin(), component_sizes.end()));

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

bool repair_and_validate_mesh(Mesh & mesh) {
	if (DEBUG) std::cout << Yellow << "      Number of removed vertices: " << ColorEnd << PMP::remove_isolated_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Number of new vertices: " << ColorEnd << PMP::duplicate_non_manifold_vertices(mesh) << std::endl;
	if (DEBUG) std::cout << Yellow << "      Number of stitches: " << ColorEnd << PMP::stitch_borders(mesh) << std::endl;
	mesh.collect_garbage();
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

void settle_mesh_z0(Mesh& mesh) {
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
}

void extrude_bottom_faces(Mesh& mesh, double target_z) {
	double z_threshold = 0.1;
	for (Vertex_index v : mesh.vertices()) {
		Point& p = mesh.point(v);
		if (p.z() >= z_threshold) {
			mesh.point(v) = Point(p.x(), p.y(), p.z() - target_z);
		}
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
	if (DEBUG) std::cout << Yellow << "      Applying translation:  " << ColorEnd << translation_vector << std::endl;
	for (auto v : mesh.vertices()) {
		mesh.point(v) = mesh.point(v) + translation_vector;
	}
}

void rotate_mesh(Mesh& mesh, double x_deg, double y_deg, double z_deg) {
	if (DEBUG) std::cout << Yellow << "      Applying Rotation:  " 
		<< ColorEnd << "(X " << x_deg << ", Y " << y_deg << ", Z " << z_deg << ")" << std::endl;
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

	
	std::cout << Yellow << "\n============================'Created by Banna'===============================" << std::endl;
	std::cout << "=============================='OCR F TOOL V3'================================\n\n" << ColorEnd << std::endl;


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

		double Width, Length, Height, Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot;
		Mesh Sub_Mesh = Result_Mesh;
		Point centroid, center;
		Result_Mesh.clear();
		
		settle_mesh_z0(Model_Mesh);
		get_dimensions(Model_Mesh, Width, Length, Height);
		//PMP::orient(Model_Mesh);
		/*get_center(Model_Mesh, center);
		translate_mesh(Model_Mesh, Kernel::Vector_3(-center.x(), -center.y() + 6, 0));*/
		get_centroid(Model_Mesh, centroid);
		translate_mesh(Model_Mesh, Kernel::Vector_3(-centroid.x(), -centroid.y() + 6, 0));

		visualize_mesh(Sub_Mesh, Model_Mesh , Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot);
		if (DEBUG) std::cout << Yellow << "      Offsets: " << ColorEnd << "X" << Model_Xoffset << ", Y" << Model_Yoffset << std::endl;
		if (DEBUG) std::cout << "               Z" << cut_height << ", Rot Z" << Model_Zrot << std::endl;

		if (Model_Xoffset != NULL || Model_Yoffset != NULL)
			translate_mesh(Model_Mesh, Kernel::Vector_3(Model_Xoffset, Model_Yoffset, 0));
		if (Model_Zrot != NULL)
			rotate_mesh(Model_Mesh, 0, 0, Model_Zrot);
		if (cut_height >= 0.1) {
			cut_mesh(Model_Mesh, cut_height, 0);
			settle_mesh_z0(Model_Mesh);
		}
		
		if (!cutting_height_Str.empty()) {
			cut_mesh(Model_Mesh, Height, std::atof(cutting_height_Str.c_str()));
			settle_mesh_z0(Model_Mesh);
		}

		if (!PMP::corefine_and_compute_union(Model_Mesh, Sub_Mesh, Result_Mesh)) {
		    std::cerr << Red << "      Model Addition failed." << ColorEnd << std::endl;
			Result_Mesh.clear();
			CGAL::copy_face_graph(Sub_Mesh, Result_Mesh);
			CGAL::copy_face_graph(Model_Mesh, Result_Mesh);
		}
	}

	if (repair_and_validate_mesh(Result_Mesh)) {
		if (DEBUG) std::cout << Green << "      Mesh repaired." << ColorEnd << std::endl;
	}
	else {
		std::cerr << Red << "      Failed to repair or validate the mesh." << ColorEnd << std::endl;
	}
	
	if (!write_STL(Output_Path_Str, Result_Mesh)) return EXIT_FAILURE;

	std::cout << Green << "      Operation completed successfully." << ColorEnd << std::endl;
	return EXIT_SUCCESS;
}