#pragma once
#ifndef BOUNDARY_VCCS_H
#define BOUNDARY_VCCS_H
#include <vector>
#include <pcl/pcl_base.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/segmentation/supervoxel_clustering.h>
#include <pcl/segmentation/lccp_segmentation.h>
#include <pcl/features/fpfh.h>

using namespace pcl;
using namespace std;

typedef PointXYZ PointT;
typedef PointCloud<PointT> PointCloudT;
typedef pair<uint32_t, int> PAIR;
typedef multimap<uint32_t, uint32_t>::iterator int_multimap;

const double eps = 1.0e-6;

class BoundaryVCCS
{
public:
	BoundaryVCCS(PointCloudT::Ptr &cloud, PointCloud<Normal>::Ptr &normal)
	{
		cloud_ = cloud;
		normalcloud_ = normal;
	}
	~BoundaryVCCS();
	//存储边界点信息
	class BoundaryData
	{
	public:
		BoundaryData(): xyz_(0.0f, 0.0f, 0.0f), normal_(0.0f, 0.0f, 0.0f)
		{
		}
		Eigen::Vector3f xyz_;
		Eigen::Vector3f normal_;
		float curvature_;
		//float distance_;
		int idx_;
	};	
	float getRvoxel()
	{
		return Rvoxel_;
	}
	void getVCCSs(float Rvoxel, float Rseed, PointCloud<PointXYZL>::Ptr &sv_labeled_cloud1, PointCloud<PointXYZL>::Ptr &sv_labeled_cloud2);
	void getBoundaryPoint(map<uint32_t, vector<BoundaryData>> &boundaries);
	void getBoundaryPoint2(map<uint32_t, vector<BoundaryData>> &boundaries);
	void expandBSV();
	void updateCentroid();
	float distCal(BoundaryData p1, Supervoxel<PointT>::Ptr p2);	
	void removeOutliers(PointCloudT::Ptr cloud1, PointCloudT::Ptr cloud2, PointCloud<Normal>::Ptr normal);
	void svGrowing_c(map < int, vector<Supervoxel<PointT>::Ptr>> &regions);
	void svGrowing_rd(map < int, vector<Supervoxel<PointT>::Ptr>> &regions);
	void seedSelection_c(vector<pair<float, int>> &curture_label);
	void seedSelection_rd(vector<pair<float, int>> &wi_label);
	int seedGrowing(int seed, int segment_number, bool flag_c);
	bool isValidated(int seed, int neighbor, bool &is_a_seed, bool flag_c);
	bool isValidated2(int seed, int neighbor, bool &is_a_seed, bool flag_c);//三种融合，效果不佳
	bool isConvex(int s_label, int t_label);
	void getRegionNeighbor();
	void mergeSmallRegions();//refinement1
	PointCloud<PointXYZL>::Ptr getLabelCloud();
	PointCloud<PointXYZL>::Ptr getPatchCloud();

	PointCloud<FPFHSignature33>::Ptr getFPFH(int label);//获取超体素的FPFH特征
private:	
	float Rvoxel_, Rseed_;
	PointCloudT::Ptr cloud_;
	PointCloud<Normal>::Ptr normalcloud_;
	map<uint32_t, Supervoxel<PointT>::Ptr > sv_clusters_;//VCCS的结果
	map<uint32_t, vector<BoundaryData>> boundaries_;//超体素标签->边界点实体
	multimap<uint32_t, uint32_t> nei_labels_;//VCCS的超体素邻接关系
	map<int, float> label_curtures_;//用于存储每个超体素的曲率均值
	map<int, float> label_wis_;//用于存储每个超体素的残差值
	map<int, int> sv_labels_;//超体素标签->区域标签
	//vector<int> num_svs_in_segment_;
	map<int, vector<Supervoxel<PointT>::Ptr>> regions_;//区域标签->内部超体素实体
	map<int, vector<int> > seglabel_to_svlist_;//区域标签->内部超体素标签
	map<int, int> svp;//存储邻接面总超体素对
	map<int, int> convex_svp;//存储邻接面凸超体素对
	map<int, set<int>> seg_label_to_neighbor_set_;//区域标签->邻接区域标签集
};

#endif BOUNDARY_VCCS_H