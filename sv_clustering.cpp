#include "sv_clustering.h"

SVclustering::~SVclustering()
{
}
void SVclustering::transFunction(PointT &p)
{
	p.x /= p.z;
	p.y /= p.z;
	p.z = log(p.z);
}

void SVclustering::selectSeedVoxel(std::vector<int> &seed_indices)
{
	pcl::octree::OctreePointCloudSearch <PointT> seed_octree(Rseed);
	seed_octree.setInputCloud(voxel_centroid_cloud_);
	seed_octree.addPointsFromInputCloud();
	std::cout << "voxel_centroid_cloud size=" << voxel_centroid_cloud_->size() << std::endl;
	std::vector<PointT, Eigen::aligned_allocator<PointT> > voxel_centers;
	int num_seeds = seed_octree.getOccupiedVoxelCenters(voxel_centers);//filtering,第一步非空鉴定
	std::cout << "Number of seed points before filtering=" << voxel_centers.size() << std::endl;
	std::vector<int> seed_indices_orig;
	seed_indices_orig.resize(num_seeds, 0);
	seed_indices.clear();//存储最终结果
	std::vector<int> closest_index;
	std::vector<float> distance;
	closest_index.resize(1, 0);
	distance.resize(1, 0);
	if (voxel_kdtree_ == 0)//初始化kd树
	{
		voxel_kdtree_.reset(new pcl::search::KdTree<PointT>);
		voxel_kdtree_->setInputCloud(voxel_centroid_cloud_);
	}
	//用KDtree搜索最靠近种子体素中心的体素
	for (int i = 0; i < num_seeds; ++i)
	{
		voxel_kdtree_->nearestKSearch(voxel_centers[i], 1, closest_index, distance);
		seed_indices_orig[i] = closest_index[0];
	}

	std::vector<int> neighbors;
	std::vector<float> sqr_distances;
	seed_indices.reserve(seed_indices_orig.size());
	float search_radius = 0.5f*Rseed;//Rsearch
	// This is 1/20th of the number of voxels which fit in a planar slice through search volume
	// Area of planar slice / area of voxel side. (Note: This is smaller than the value mentioned in the original paper)
	float min_points = 0.05f * (search_radius)*(search_radius)* 3.1415926536f / (Rvoxel*Rvoxel);
	for (size_t i = 0; i < seed_indices_orig.size(); ++i)
	{
		int flag;
		int num = voxel_kdtree_->radiusSearch(seed_indices_orig[i], search_radius, neighbors, sqr_distances);
		int min_index = seed_indices_orig[i];
		if (num > min_points)
		{
			//对非孤立点再判断是否为边界
			flag = isboundaryIndex(min_index);
			if (flag)//对非边界点体素索引保存
			{
				seed_indices.push_back(min_index);
			}
			else//否则对边界点邻域索引判断
			{
				// flag = 1;
				for (int j = 1; j < neighbors.size(); j++)
				{
					flag = isboundaryIndex(neighbors[j]);
					if (flag)
					{
						seed_indices.push_back(neighbors[j]);
						break;
					}
				}
			}
		}

	}
	std::cout << "seed_indices size=" << seed_indices.size() << std::endl;
}
//voxelization
//[in] cloud&boundary
void SVclustering::computeSupervoxel(const PointCloudT::ConstPtr &cloud, const CloudNormal::ConstPtr &normal_cloud, CloudBoundary &boundary)
{
	/*Rvoxel = 0.007f;
	Rseed = 0.055f;*/
	float Wcolor = 0.0f;
	float Wspatial = 2.5f;
	float Wnormal = 5.5f;
	//bool use_single_cam_transform = true;//logZ
	//pcl::SupervoxelClustering<PointT> super(Rvoxel, Rseed);
	//std::map<uint32_t, pcl::Supervoxel<PointT>::Ptr> supervoxel_clusters;
	//multimap<uint32_t, uint32_t> supervoxel_adjacency;
	//super.setUseSingleCameraTransform(use_single_cam_transform);
	//super.setInputCloud(cloud);
	//super.setNormalCloud(normals);
	adjacency_octree.reset(new OctreeAdjacencyT(Rvoxel));//要先初始化，不然共享指针会报错
	adjacency_octree->setInputCloud(cloud);
	//adjacency_octree->setTransformFunction(boost::bind(&transFunction, adjacency_octree, _1));
	adjacency_octree->addPointsFromInputCloud();
	//std::cout << "Size of octree =" << adjacency_octree->getLeafCount() << std::endl;

	voxel_centroid_cloud_.reset(new PointCloudT);
	voxel_centroid_cloud_->resize(adjacency_octree->getLeafCount());
	std::vector<PointT, Eigen::aligned_allocator<PointT> > voxel_center_list_arg;
	adjacency_octree->getOccupiedVoxelCenters(voxel_center_list_arg);
	LeafVectorT::iterator leaf_itr = adjacency_octree->begin();
	//PointCloudT::iterator cent_cloud_itr = voxel_centroid_cloud_->begin();
	for (int idx = 0; leaf_itr != adjacency_octree->end(); ++leaf_itr, ++idx)
	{
		VoxelData& new_voxel_data1 = (*leaf_itr)->getData();//VoxelData未赋值
		voxel_centroid_cloud_->points[idx] = voxel_center_list_arg[idx];
		new_voxel_data1.xyz_ = voxel_center_list_arg[idx].getVector3fMap();
		new_voxel_data1.idx_ = idx;
	}

	
	for (int i = 0; i<cloud->size(); i++)
	{
		if (!pcl::isFinite<PointT>(cloud->points[i]))
			continue;
		LeafContainerT* leaf;
		leaf = adjacency_octree->getLeafContainerAtPoint(cloud->points[i]);
		VoxelData& new_voxel_data = leaf->getData();
		new_voxel_data.normal_ += normal_cloud->at(i).getNormalVector4fMap();
		new_voxel_data.curvature_ += normal_cloud->at(i).curvature;
		if (boundary[i].boundary_point>0)
		{
			int index = new_voxel_data.idx_;
			bounday_index[i] = index;
		}

	}
	for (leaf_itr = adjacency_octree->begin(); leaf_itr != adjacency_octree->end(); ++leaf_itr)
	{
		VoxelData& new_voxel_data = (*leaf_itr)->getData();
		new_voxel_data.normal_.normalize();
		new_voxel_data.owner_ = 0;
		new_voxel_data.distance_ = std::numeric_limits<float>::max();//编译器允许float类型数最大值
		int num_points = (*leaf_itr)->getPointCounter();
		new_voxel_data.curvature_ /= num_points;
	}

	/*super.setColorImportance(Wcolor);
	super.setSpatialImportance(Wspatial);
	super.setNormalImportance(Wnormal);
	super.extract(supervoxel_clusters);
	std::cout << "Found " << supervoxel_clusters.size() << " Supervoxels!\n";*/
	//PointCloud<PointXYZL>::Ptr labeled_voxel_cloud = super.getLabeledVoxelCloud();
	//showCloud(labeled_voxel_cloud);
	//super.getSupervoxelAdjacency(supervoxel_adjacency);
}
int SVclustering::isboundaryIndex(int mindex)
{
	std::map<int, int>::iterator iter;
	for (iter = bounday_index.begin(); iter != bounday_index.end(); iter++)
	{
		if (mindex == iter->second)//为边界点体素，需要继续判断相邻体素
		{
			return 0;
			//std::cout << "it is boundary seed!!!" << std::endl;
		}

	}
	return 1;
}
void SVclustering::createSVHelpers(std::vector<int> &seed_indices)
{
	supervoxel_helpers_.clear();
	for (size_t i = 0; i < seed_indices.size(); ++i)
	{
		//label_+parent_初始化
		supervoxel_helpers_.push_back(new SVHelper(i + 1, this));
		//Find which leaf corresponds to this seed index
		LeafContainerT* seed_leaf = adjacency_octree->at(seed_indices[i]);//adjacency_octree_->getLeafContainerAtPoint (seed_points[i]);
		if (seed_leaf)
		{
			supervoxel_helpers_.back().addLeaf(seed_leaf);//添加seed_leaf指针,ower_=this
		}
		else
		{
			PCL_WARN("Could not find leaf in pcl::SupervoxelClustering<PointT>::createSupervoxelHelpers - supervoxel will be deleted \n");
		}
	}
}

void SVclustering::SVHelper::addLeaf(LeafContainerT* leaf_arg)
{
	leaves_.insert(leaf_arg);//LeafSetT leaves_;
	voxeldata& voxel_data = leaf_arg->getData();
	//voxel_data.owner_ = this;
}
