#pragma once
#ifndef SV_CLUSTERING_H
#define SV_CLUSTERING_H
#include "pre_boundary.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/segmentation/supervoxel_clustering.h>
#include <pcl/octree/octree_search.h>
#include <pcl/octree/octree_pointcloud.h>
#include <pcl/octree/octree_pointcloud_adjacency.h>

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
typedef pcl::PointCloud<pcl::Normal> CloudNormal;
typedef pcl::PointCloud<pcl::Boundary> CloudBoundary;
//typedef pcl::SupervoxelClustering<PointT>::OctreeAdjacencyT OctreeAdjacencyTT;
//typedef pcl::SupervoxelClustering<PointT>::VoxelData VoxelData;
//typedef pcl::octree::OctreePointCloudAdjacencyContainer<PointT, VoxelData> LeafContainerT;
//typedef std::vector <LeafContainerT*> LeafVectorT;
//typedef pcl::octree::OctreePointCloudAdjacency<PointT, LeafContainerT> OctreeAdjacencyT;
//typedef pcl::octree::OctreePointCloudSearch <PointT> OctreeSearchT;

class SVclustering : public pcl::SupervoxelClustering<PointT>
{
public:
	SVclustering(float voxel, float seed) : pcl::SupervoxelClustering<PointT>(voxel, seed)
	{
		Rvoxel = voxel;
		Rseed = seed;
	};
	~SVclustering();
	void transFunction(PointT &p);
	//void computeLeafIdx(SupervoxelClustering<PointT>::OctreeAdjacencyT::Ptr &adjacency_octree);
	void computeSupervoxel(const PointCloudT::ConstPtr &cloud, const CloudNormal::ConstPtr &normal_cloud, CloudBoundary &boundary);
	void selectSeedVoxel(std::vector<int> &seed_indices);
	int isboundaryIndex(int min_index);
	void createSVHelpers(std::vector<int> &seed_indices);

	class SVHelper;
	friend class SVHelper;
	class voxeldata : public VoxelData
	{
		SVHelper* owner_;
	};
	typedef pcl::octree::OctreePointCloudAdjacencyContainer<PointT, voxeldata> LeafContainerT;
	typedef std::vector <LeafContainerT*> LeafVectorT;
	typedef pcl::octree::OctreePointCloudAdjacency<PointT, LeafContainerT> OctreeAdjacencyT;
	OctreeAdjacencyT::Ptr adjacency_octree;
	pcl::search::KdTree<PointT>::Ptr voxel_kdtree_;
	PointCloudT::Ptr voxel_centroid_cloud_;
	std::map<int, int> bounday_index;
	
	class SVHelper
	{
	public:
		struct compareLeaves
		{
			bool operator() (LeafContainerT* const &left, LeafContainerT* const &right) const
			{
				const VoxelData& leaf_data_left = left->getData();
				const VoxelData& leaf_data_right = right->getData();
				return leaf_data_left.idx_ < leaf_data_right.idx_;
			}
		};
		typedef std::set<LeafContainerT*, SVHelper::compareLeaves> LeafSetT;//×Ô¶¨ÒåÅÅÐòsetÈÝÆ÷
		typedef LeafSetT::iterator iterator;
		typedef LeafSetT::const_iterator const_iterator;
		SVHelper(uint32_t label, SupervoxelClustering* parent_arg) :
			label_(label),
			parent_(parent_arg)
		{ }
		void addLeaf(LeafContainerT* leaf_arg);
	private:
		//Stores leaves
		LeafSetT leaves_;
		uint32_t label_;
		VoxelData centroid_;
		SupervoxelClustering* parent_;
	};
	
	friend void boost::checked_delete<>(const SVHelper *);
	typedef boost::ptr_list<SVHelper> HelperListT;
	HelperListT supervoxel_helpers_;
	
private:
	float Rvoxel;
	float Rseed;

};

#endif SV_CLUSTERING_H /*SV_CLUSTERING_H*/
