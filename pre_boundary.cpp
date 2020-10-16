#include "pre_boundary.h"

PreBoundary::~PreBoundary()
{
}

void PreBoundary::computeNormals(const PointCloudT::ConstPtr &cloud,float r=0.02)
{
	pcl::NormalEstimation<PointT, pcl::Normal> normalEst;
	CloudNormal::Ptr normals(new CloudNormal());
	normalEst.setInputCloud(cloud);
	//normalEst.setKSearch(100);
	normalEst.setRadiusSearch(r);
	normalEst.compute(*normals);
	normals_ = normals;
}
//利用PCL自带的边界计算
PointCloudT::Ptr PreBoundary::computeCloudBoundary(const PointCloudT::ConstPtr &cloud)
{
	pcl::BoundaryEstimation<PointT, pcl::Normal, pcl::Boundary> boundaryEst;
	PointCloudT::Ptr cloud_boundary(new PointCloudT());
	CloudBoundary boundaries;
	boundaryEst.setInputCloud(cloud);
	boundaryEst.setInputNormals(normals_);
	//boundaryEst.setKSearch(200);
	boundaryEst.setRadiusSearch(0.02);
	boundaryEst.setAngleThreshold(M_PI / 2);
	boundaryEst.setSearchMethod(pcl::search::KdTree<PointT>::Ptr(new pcl::search::KdTree<PointT>()));
	boundaryEst.compute(boundaries);
	boundaries_ = boundaries;
	//save boundaries as PointXYZRGBA type
	for (int i = 0; i < cloud->points.size(); i++)
	{
		if (boundaries[i].boundary_point>0)
		{
			cloud_boundary->push_back(cloud->points[i]);
		}
	}
	return cloud_boundary;
}
//利用均值漂移和局部密度获取边界
PointCloudT::Ptr PreBoundary::computeCloudBoundary2(const PointCloudT::ConstPtr &cloud, int k, double m)
{
	PointCloudT::Ptr cloud_boundary(new PointCloudT());
	double min_dis = 0.0;
	std::vector<int> indices(k);
	std::vector<float> sqr_distances(k);
	PointT centroid;
	PointT q_point;
	pcl::search::KdTree<PointT> tree;
	//pcl::KdTreeFLANN<PointT> tree;
	tree.setInputCloud(cloud);
	for (size_t i = 0; i < cloud->size(); ++i)
	{
		q_point = cloud->points[i];
		/*if (!isfinite((*cloud)[i].x))
		{
			continue;
		}*/
		if (tree.nearestKSearch(i, k, indices, sqr_distances) == k)
		{
			min_dis = sqr_distances[1];
			for (int i = 2; i < k; i++)
			{
				if (min_dis>sqr_distances[i])
				{
					min_dis = sqr_distances[i];
				}
			}
			for (int i = 0; i < k; i++)
			{
				centroid.x += cloud->points[indices[i]].x;
				centroid.y += cloud->points[indices[i]].y;
				centroid.z += cloud->points[indices[i]].z;
			}
			centroid.x /= k;
			centroid.y /= k;
			centroid.z /= k;
			if (((centroid.x - q_point.x)*(centroid.x - q_point.x) + (centroid.y - q_point.y)*(centroid.y - q_point.y) + (centroid.z - q_point.z)*(centroid.z - q_point.z))>m*m*min_dis)
			{
				cloud_boundary->push_back(q_point);
			}

		}

	}
	return cloud_boundary;
}
//论文2685的部分复现，缺少特征线tracing，效果一般
PointCloudT::Ptr PreBoundary::computeCloudBoundaryAGPN(const PointCloudT::ConstPtr &cloud, int k, double distance)
{
	pcl::BoundaryEstimation<PointT, pcl::Normal, pcl::Boundary> boundaryEst;
	CloudNormal::Ptr normals(new CloudNormal());
	pcl::PointCloud<pcl::Boundary> boundaries;
	boundaries.resize(cloud->size());
	PointCloudT::Ptr cloud_boundary(new PointCloudT());
	pcl::search::KdTree<PointT> tree;
	tree.setInputCloud(cloud);
	std::vector<int> indices(k);
	std::vector<float> sqr_distances(k);
	Eigen::VectorXf model_coefficients;

	for (int i = 0; i < cloud->size(); i++)
	{
		pcl::Normal p;
		/*if (!isfinite((*cloud)[i].x))
		{
			p.normal_x = p.normal_y = p.normal_z = p.curvature = std::numeric_limits<float>::quiet_NaN();
			normals->push_back(p);
			boundaries[i].boundary_point = std::numeric_limits<uint8_t>::quiet_NaN();
			continue;
		}*/
		if (tree.nearestKSearch(i, k, indices, sqr_distances) == k)
		{
			pcl::SampleConsensusModelPlane<PointT>::Ptr model_plane(new pcl::SampleConsensusModelPlane<PointT>(cloud, indices));
			std::vector<int> inliers;
			pcl::RandomSampleConsensus<PointT> ransac(model_plane);
			ransac.setDistanceThreshold(distance);
			ransac.computeModel();
			ransac.getModelCoefficients(model_coefficients);
			ransac.getInliers(inliers);
	
			p.normal_x = model_coefficients[0];
			p.normal_y = model_coefficients[1];
			p.normal_z = model_coefficients[2];
			p.curvature = 0;
			normals->push_back(p);
			std::vector<int>::iterator result = find(inliers.begin(), inliers.end(), i);
			if (result == inliers.end() || inliers.size() < 3)
			{
				boundaries[i].boundary_point = std::numeric_limits<uint8_t>::quiet_NaN();
				continue;
			}
			PointT q_point = cloud->points[i];
			Eigen::Vector4f u = Eigen::Vector4f::Zero(), v = Eigen::Vector4f::Zero();
			pcl::Vector4fMap p_coeff_v = p.getNormalVector4fMap();
			v = p_coeff_v.unitOrthogonal();
			u = p_coeff_v.cross3(v);
			float angle_threshold = M_PI / 2;
			boundaries[i].boundary_point = isBoundary(cloud, q_point, inliers, u, v, angle_threshold);
			continue;	
		}
		boundaries[i].boundary_point = std::numeric_limits<uint8_t>::quiet_NaN();
		p.normal_x = p.normal_y = p.normal_z = p.curvature = std::numeric_limits<float>::quiet_NaN();
		normals->push_back(p);
	}
	normals_ = normals;
	boundaries_ = boundaries;
	for (int i = 0; i < cloud->points.size(); i++)
	{
		if (boundaries[i].boundary_point>0)
		{
			cloud_boundary->push_back(cloud->points[i]);
		}
	}
	//pcl::StatisticalOutlierRemoval<PointT> Static;   //创建滤波器对象
	//Static.setInputCloud(cloud_boundary);                           //设置待滤波的点云
	//Static.setMeanK(15);                               //设置在进行统计时考虑查询点临近点数
	//Static.setStddevMulThresh(0.05);                      //设置判断是否为离群点的阀值
	//Static.filter(*cloud_boundary);                    //存储
	return cloud_boundary;
}
bool PreBoundary::isBoundary(const PointCloudT::ConstPtr &cloud, PointT q_point, std::vector<int> &indices, Eigen::Vector4f &u, const Eigen::Vector4f &v, float angle_threshold)
{
	if (indices.size() < 3)
		return (false);
	/*if (!pcl_isfinite(q_point.x) || !pcl_isfinite(q_point.y) || !pcl_isfinite(q_point.z))
		return (false);*/
	std::vector<float> angles(indices.size());
	float max_dif = FLT_MIN, dif;
	int cp = 0;
	for (size_t i = 0; i < indices.size(); ++i)
	{
		/*if (!pcl_isfinite(cloud->points[indices[i]].x) || !pcl_isfinite(cloud->points[indices[i]].y) || !pcl_isfinite(cloud->points[indices[i]].z))
			continue;*/

		Eigen::Vector4f delta = cloud->points[indices[i]].getVector4fMap() - q_point.getVector4fMap();
		if (delta == Eigen::Vector4f::Zero())
			continue;

		angles[cp++] = atan2f(v.dot(delta), u.dot(delta)); // the angles are fine between -PI and PI too与X轴方向角
	}
	if (cp == 0)
		return (false);

	angles.resize(cp);
	std::sort(angles.begin(), angles.end());
	//find max angel value
	for (size_t i = 0; i < angles.size() - 1; ++i)
	{
		dif = angles[i + 1] - angles[i];
		if (max_dif < dif)
			max_dif = dif;
	}
	dif = 2 * static_cast<float> (M_PI)-angles[angles.size() - 1] + angles[0];
	if (max_dif < dif)
		max_dif = dif;

	// Check results
	if (max_dif > angle_threshold)
		return (true);
	else
		return (false);

}