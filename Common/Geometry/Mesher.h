#pragma once
#include <vector>
#include <string>
#include <memory>
#include "GeometryPrimitives.h"

namespace FemLibrary {

	struct SmoothingParameters
	{
		bool smoothing = true;
		double smoothingIterationsMultiplier = 5; // 0 means no smoothing
		double minQualitySmoothingRatio = 0.1; // 0 <= value <= 1, 0.1 means 10% - minQuality, 90% - meanQuality
		double moveThreshold = 0.00001; // minimal acceptable move length = moveThreshold * ribLength(oldPoint)
	};

	struct MeshingResult
	{
		bool success;
		std::vector<Vector2D> points;
		std::vector<int> triangleIndexes;
		std::string errorMsg;
		void Clear() { points.clear(); triangleIndexes.clear(); errorMsg.clear(); success = false; }
	};

	struct MeshingParameters
	{
		SmoothingParameters smoothing_params;

		int number_of_integration_steps = 100;
		int boundary_optimization_iterations = 20;

		bool ConvexDelaynayOnly = false; // for unittests
								 
		double ToleranceRatio = 1.5;
		int MaxNumberOfTriangles = 0;

		//
		double RibLength = 1; // replace with function

		// intermediate states
		bool saveIntermediateSteps = false;

		bool extendedStatistics = false;
	};
	
	struct MeshingStatistics {
		static const int displaySize = 200;
		bool isExtended = false;
		
		__int64 startTime = 0, endTime = 0;
		double timeInputProcessing    = .0;
		double timeConvexDelaunay     = .0;
		double timeConstainedDelaunay = .0;
		double timeAFM                = .0;
		double timeSmoothingMin       = .0;
		double timeSmoothingMean      = .0;
		double timeTotal              = .0;

		int numPoints    = 0;
		int numTriangles = 0;
		
		double minAngle   = .0, maxAngle   = .0;
		double minQuality = .0, maxQuality = .0;
		std::vector<Vector2D> histTriangleAngles;
		std::vector<Vector2D> histTriangleQualities;
		std::vector<Vector2D> histTriangleAnglesAO1;
		std::vector<Vector2D> histTriangleQualitiesAO1;
		std::vector<Vector2D> histTriangleAnglesAO2;
		std::vector<Vector2D> histTriangleQualitiesAO2;
					
		std::vector<Vector2D>  minQualityDuringOptimization;
		std::vector<Vector2D> meanQualityDuringOptimization;
		int minQualityIterations = 0;
		int meanQualityIterations = 0;

		int  minQualityMoveFailureCounter = 0;
		int meanQualityMoveFailureCounter = 0;

		int  minQualityTooSmalllMoveCounter = 0;
		int meanQualityTooSmalllMoveCounter = 0;

		int  minQualitySuccessfulMoveCounter = 0;
		int meanQualitySuccessfulMoveCounter = 0;
	};

	void SaveMeshingStatisticsToTxtFile(const MeshingStatistics& stat, const char* filename);

	struct IntermediateStep
	{
		Vector2D newPoint;
		std::vector<SimpleTriangle> triangles;
	};

	class FEMLIBRARY_API Mesh
	{
		class MeshImpl;
		std::unique_ptr<MeshImpl> d;

	public:
		Mesh();
		~Mesh();
		Mesh(Mesh const&) = delete;
		Mesh& operator=(Mesh const&) = delete;


		MeshingResult generateMesh(const FragmentGeometry& cg, const MeshingParameters& params) const;
		std::vector<IntermediateStep> getIntermediateSteps() const;
		MeshingStatistics getStatistics() const;
	};
};