#include "boundary_VCCS.h"
#include "common.h"

BoundaryVCCS::~BoundaryVCCS()
{
}
//��ȡ�����أ�δʹ����ɫ��Ϣ
void BoundaryVCCS::getVCCSs(float Rvoxel, float Rseed, PointCloud<PointXYZL>::Ptr &sv_labeled_cloud1, PointCloud<PointXYZL>::Ptr &sv_labeled_cloud2)
{
	Rvoxel_ = Rvoxel;
	Rseed_ = Rseed;
	SupervoxelClustering<PointT> super(Rvoxel, Rseed);
	super.setInputCloud(cloud_);
	//Ȩ��ΪĬ��ֵ���ɸ���ʵ�ʳ�������
	float wc = 0.0f;
	float ws = 1.0f;//ws=2;wn=1
	float wn = 4.0f;//Ĭ��ws=1,wn=4
	super.setColorImportance(wc);
	super.setSpatialImportance(ws);
	super.setNormalImportance(wn);
	super.setNormalCloud(normalcloud_);
	//super.setUseSingleCameraTransform(true);//omit
	super.extract(sv_clusters_);
	/*for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end();)
	{
		PointCloudT::Ptr voxels = labelsv_itr->second->voxels_;
		if (voxels->size() == 1)
		{
			sv_clusters_.erase(labelsv_itr++);
		}
		else
		{
			labelsv_itr++;
		}
	}*/
	super.getSupervoxelAdjacency(nei_labels_);//get
	sv_labeled_cloud1 = super.getLabeledVoxelCloud();
	//����
	//int max_depth = static_cast<int> (1.8f*Rseed / Rvoxel);
	for (int i = 0; i < 10; i++)
	{
		map<uint32_t, vector<BoundaryVCCS::BoundaryData>> cloud_boundary;
		getBoundaryPoint2(cloud_boundary);
		expandBSV();
		//����sv_clusters_������
		updateCentroid();
	}
	sv_labeled_cloud2 = getLabelCloud();
}
//1.��ֱͶӰ�����õĵ�������΢�ٵ�
void BoundaryVCCS::getBoundaryPoint(map<uint32_t, vector<BoundaryData>> &boundaries)
{
	//����supervoxel_clusters_(lable_sv)
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		uint32_t label = labelsv_itr->first;
		Supervoxel<PointT>::Ptr cur_sv = labelsv_itr->second;
		Eigen::Vector3f cur_centroid = cur_sv->centroid_.getVector3fMap();
		Eigen::Vector3f cur_normal = cur_sv->normal_.getNormalVector3fMap().normalized();
		PointCloudT::Ptr voxels = cur_sv->voxels_;
		PointCloud<Normal>::Ptr normals = cur_sv->normals_;
		vector<BoundaryData> pns;
		for (int i = 0; i < voxels->size(); i++)
		{
			Eigen::Vector3f q_point = voxels->points[i].getVector3fMap();
			Eigen::Vector3f delta = q_point - cur_centroid;
			if (delta == Eigen::Vector3f::Zero())
				continue;
			Eigen::Vector3f q_proj = delta - cur_normal*(delta.dot(cur_normal));
			float Dproj = q_proj.dot(q_proj);//����ƽ����
			if (sqrt(Dproj) > Rseed_ / 2)
			{
				BoundaryData pn;
				pn.xyz_ = voxels->points[i].getVector3fMap();
				pn.normal_ = normals->points[i].getNormalVector3fMap();
				pn.idx_ = i;
				pns.push_back(pn);
			}
			else
				continue;
		}
		boundaries[label] = pns;
	}
	boundaries_ = boundaries;
	//cout << "boundaries size = " << boundaries->size() << endl;
}
//2.������Ʋ�ؾ��룬͸��ͶӰ�������,��������
void BoundaryVCCS::getBoundaryPoint2(map<uint32_t, vector<BoundaryData>> &boundaries)
{
	//map<uint32_t, float> curvatures;
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		uint32_t label = labelsv_itr->first;
		Supervoxel<PointT>::Ptr cur_sv = labelsv_itr->second;
		Eigen::Vector3f pc = cur_sv->centroid_.getVector3fMap();
		PointCloudT::Ptr voxels = cur_sv->voxels_;
		PointCloud<Normal>::Ptr normals = cur_sv->normals_;
		vector<BoundaryData> pns;
		//PointCloudT::Ptr boundaryCloud(new PointCloudT);
		float curvature = 0;
		for (int i = 0; i < normals->size(); i++)
		{
			curvature += normals->points[i].curvature;
		}
		//curvatures[label] = sum_curvature / cur_normals->size();//��ȡÿ��SV�����ʣ�ƽ����
		curvature = curvature / normals->size();
		float r = 1 / curvature;//���ʰ뾶
		for (int i = 0; i < voxels->size(); i++)
		{
			Eigen::Vector3f pi = voxels->points[i].getVector3fMap();//������Ϊԭ��
			pi = pi - pc;//ȥ�Ļ�
			float Dproj = r*sqrt(pi[0] * pi[0] + pi[1] * pi[1] + 4 * pi[2] * pi[2]) / abs(r - pi[2]);
			if (Dproj > Rseed_ / 2)
			{
				BoundaryData pn;
				pn.xyz_ = voxels->points[i].getVector3fMap();
				pn.normal_ = normals->points[i].getNormalVector3fMap();
				pn.curvature_ = normals->points[i].curvature;
				pn.idx_ = i;
				pns.push_back(pn);
			}
			else
				continue;
		}
		boundaries[label] = pns;
	}
	boundaries_ = boundaries;
}
//�߽��������SV���ľ��룬�����׼���Ը���������
float BoundaryVCCS::distCal(BoundaryData p1, Supervoxel<PointT>::Ptr p2)
{
	float wn = 1;
	float wd = 0.1/ Rvoxel_;//�������Ȩ�ؽϴ��ʹ�����ؿ�Խ�߽�
	Eigen::Vector3f dnormal = p1.normal_ - p2->normal_.getNormalVector3fMap();
	Eigen::Vector3f dxyz = p1.xyz_ - p2->centroid_.getVector3fMap();
	return wn*dnormal.norm() + wd*dxyz.norm();
}
//K��ֵ���ࣺ�߽��ͳ��������ĵ�
void BoundaryVCCS::expandBSV()
{
	map<uint32_t, vector<BoundaryData>>::iterator b_itr;
	for (b_itr = boundaries_.begin(); b_itr != boundaries_.end(); ++b_itr)
	{
		uint32_t label = b_itr->first;//��ǰ�����ر�ǩ
		PointCloud<PointT>::Ptr outliers(new PointCloud<PointT>);
		pair<int_multimap, int_multimap> cur_labels = nei_labels_.equal_range(label);//����sv
		vector<BoundaryData> curBoundaries = b_itr->second;
		Supervoxel<PointT>::Ptr pc0(new Supervoxel<PointT>);
		pc0 = sv_clusters_[label];//��ǰ�߽��ĳ�����
		//ѭ����ǰ�����������б߽��
		for (int i = 0; i < curBoundaries.size(); i++)
		{
			//����ÿ���߽������Χ�����������ĵ����С���룬����¼��Ӧ�ı�ǩ
			BoundaryData pn = curBoundaries[i];	//��ǰ�߽������
			uint32_t min_label = label;
			float min_dist = distCal(pn, pc0);//�뵱ǰ�����ؾ�����Ϊ��Сֵ
			for (auto iter = cur_labels.first; iter != cur_labels.second; ++iter)
			{
				uint32_t cur_label = iter->second;
				Supervoxel<PointT>::Ptr pc = sv_clusters_[cur_label];
				if (distCal(pn, pc) < min_dist)
				{
					min_dist = distCal(pn, pc);
					min_label = cur_label;//�ҵ�������С����ĳ�����
				}
			}
			//���ڲ������ʺϵı߽�����
			if (min_label != label)
			{			
				PointT p;
				p.x = pn.xyz_[0];
				p.y = pn.xyz_[1];
				p.z = pn.xyz_[2];
				Normal n;
				n.normal_x = pn.normal_[0];
				n.normal_y = pn.normal_[1];
				n.normal_z = pn.normal_[2];
				n.curvature = pn.curvature_;
				sv_clusters_[min_label]->voxels_->push_back(p);	
				sv_clusters_[min_label]->normals_->push_back(n);
				outliers->push_back(p);//�洢��ǰsvҪɾ���ı߽��
			}
		}
		//�Ե�ǰ�������ѵ������ı߽�㣨����&���ߣ�����ɾ��
		removeOutliers(outliers, pc0->voxels_, pc0->normals_);		
	}

}
//update sv_clusters_ information for the new interation
void BoundaryVCCS::updateCentroid()
{
	float cenDist = 0;
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		Supervoxel<PointT>::Ptr clusters = labelsv_itr->second;
		PointCloudT::Ptr voxels = clusters->voxels_;
		PointCloud<Normal>::Ptr normals = clusters->normals_;
		Eigen::Vector3f centroid=Eigen::Vector3f::Zero();
		Eigen::Vector3f normal = Eigen::Vector3f::Zero();
		for (int i = 0; i < voxels->size(); i++)
		{
			normal += normals->points[i].getNormalVector3fMap();
			centroid += voxels->points[i].getVector3fMap();
		}
		normal.normalize();
		centroid /= voxels->size();
		Eigen::Vector3f pre_centroid = clusters->centroid_.getVector3fMap();
		cenDist += (centroid - pre_centroid).norm();
		clusters->centroid_.x = centroid[0];
		clusters->centroid_.y = centroid[1];
		clusters->centroid_.z = centroid[2];
		clusters->normal_.normal_x = normal[0];
		clusters->normal_.normal_y = normal[1];
		clusters->normal_.normal_z = normal[2];
	}
	cout << "����ǰ�����Ĳ�ƽ��֮�ͣ�" << cenDist << endl;
}
//ɾ��cloud2�к���cloud1�ĵ�
void BoundaryVCCS::removeOutliers(PointCloudT::Ptr cloud1, PointCloudT::Ptr cloud2, PointCloud<Normal>::Ptr normal)
{
	KdTreeFLANN<PointT> kdtree;
	PointT searchPoint;
	int K = 1;
	vector<int> pointIdxNKNSearch(K);      //�洢��ѯ���������
	vector<float> pointNKNSquaredDistance(K); //�洢���ڵ��Ӧ����ƽ��
	int num = 0;
	//vector<PointT> DeleteData;
	for (auto iter1 = cloud1->begin(); iter1 != cloud1->end(); ++iter1)
	{
		searchPoint.x = iter1->x;
		searchPoint.y = iter1->y;
		searchPoint.z = iter1->z;
		kdtree.setInputCloud(cloud2);
		num = kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);//��cloud2���ҵ�cloud1�ж�Ӧ�ĵ㣨��Ӧcloud2��index��
		if (num > 0)
		{
			if (sqrt(pointNKNSquaredDistance[0]) < eps)
			{
				auto iter2 = cloud2->begin() + pointIdxNKNSearch[0];
				auto iter3 = normal->begin() + pointIdxNKNSearch[0];
				cloud2->erase(iter2);//delete point xyz_
				normal->erase(iter3);//delete point normal_
				//DeleteData.push_back(searchPoint);
				if (cloud2->size() == 0)
				{
					break;
				}
				//reset
				searchPoint.x = 0;
				searchPoint.y = 0;
				searchPoint.z = 0;
				num = 0;
				pointIdxNKNSearch.clear();
				pointNKNSquaredDistance.clear();
			}
		}
	}
	//cout << DeleteData.size() << endl;
}
//��������������1:������ѡȡ���ӵ㣬smoothness��Ϊ����׼�򣬻��facet
void BoundaryVCCS::seedSelection_c(vector<pair<float, int>> &curture_label)
{
	//int i_point = 0;
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		int label = labelsv_itr->first;
		sv_labels_[label] = 0;//��ʼ��
		PointCloud<Normal>::Ptr normals = labelsv_itr->second->normals_;
		float curvature = 0;
		int count = normals->size();
		for (int i = 0; i < count; i++)
		{
			float tmp = normals->points[i].curvature;
			if (isnan(tmp))
			{
				count--;
				continue;
			}
			curvature += tmp;
		}
		curvature /= count;
		label_curtures_[label] = curvature;
		if (isnan(curvature))
		{
			continue;
		}
		curture_label.push_back(pair<float, int>(curvature, label));
	}
	sort(curture_label.begin(), curture_label.end());
}
void BoundaryVCCS::svGrowing_c(map<int, vector<Supervoxel<PointT>::Ptr>> &regions)
{
	int num_of_sv = sv_clusters_.size();
	//����ÿ������������ֵ������������
	vector<pair<float, int>> curture_label;//�������򣨰���ֵ�����Ϊ���㣩
	seedSelection_c(curture_label);
	//��������ѡ�����ӽ�������
	int seed_counter = 0;//get next seed
	int seed = curture_label[seed_counter].second;//minmum curture -> label
	int segmented_num = 0;//�Ѿ��ָ��sv��
	int segment_label = 1;//�ָ��ǩ
	while (segmented_num < num_of_sv)
	{
		int sv_in_segment;
		sv_in_segment = seedGrowing(seed, segment_label, true);
		segmented_num += sv_in_segment;
		//num_svs_in_segment_.push_back(sv_in_segment);
		segment_label++;
		//find the next seed that is not segmented
		for (int i_seed = seed_counter + 1; i_seed < curture_label.size(); i_seed++)
		{
			int index = curture_label[i_seed].second;
			if (sv_labels_[index] == 0)
			{
				seed = index;
				seed_counter = i_seed;
				break;
			}
		}
	}
	//����ǩһ����sv�ϲ�
	regions.clear();
	for (map<int, int>::iterator itr = sv_labels_.begin(); itr != sv_labels_.end(); ++itr)
	{
		int label = itr->second;
		int svlabel = itr->first;
		if (label == 0)
			continue;
		regions[label].push_back(sv_clusters_[svlabel]);	
		seglabel_to_svlist_[label].push_back(svlabel);
	}
	regions_ = regions;
}
//��������������2���Բв�ֵ���ܶ�ѡȡ���ӵ㣬����׼�����
void BoundaryVCCS::seedSelection_rd(vector<pair<float, int>> &wi_label)
{
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		int label = labelsv_itr->first;
		sv_labels_[label] = 0;//��ʼ��
		PointCloudT::Ptr voxels = labelsv_itr->second->voxels_;
		float residual = 0;
		int num = voxels->size();
		float density;
		Eigen::Vector3f pc = labelsv_itr->second->centroid_.getVector3fMap();
		Eigen::Vector3f nc = labelsv_itr->second->normal_.getNormalVector3fMap();
		for (int i = 0; i < num; i++)
		{
			Eigen::Vector3f pi = voxels->at(i).getVector3fMap();
			residual += std::pow(nc.dot(pi - pc), 2);
		}
		residual = sqrt(residual / voxels->size());
		//label_residuals_[label] = residual;
		density = 1.0 / num;
		float k = 0.1;
		//float wi = 1 - exp(-(residual*residual) / (2 * k*k)) * exp(-(density*density) / (2 * k*k));
		float wi = residual;
		label_wis_[label] = wi;
		if (isnan(wi))
		{
			continue;
		}
		wi_label.push_back(pair<float, int>(wi, label));
	}
	sort(wi_label.begin(), wi_label.end());

}
void BoundaryVCCS::svGrowing_rd(map < int, vector<Supervoxel<PointT>::Ptr>> &regions)
{
	//���㳬�����ڲ��в�ֵ
	int num_of_sv = sv_clusters_.size();
	vector<pair<float, int>> wi_label;
	seedSelection_rd(wi_label);
	
	//����wiѡȡ���ӣ��в�ֵ������С�����������ܶ�
	int seed_counter = 0;//get next seed
	int seed = wi_label[seed_counter].second;//minmum rediual+1/num -> label
	int segmented_num = 0;//�Ѿ��ָ��sv��
	int segment_label = 1;//�ָ��ǩ
	while (segmented_num < num_of_sv)
	{
		int sv_in_segment;
		sv_in_segment = seedGrowing(seed, segment_label, false);
		segmented_num += sv_in_segment;
		//num_svs_in_segment_.push_back(sv_in_segment);
		segment_label++;
		//find the next seed that is not segmented
		for (int i_seed = seed_counter + 1; i_seed < wi_label.size(); i_seed++)
		{
			int index = wi_label[i_seed].second;
			if (sv_labels_[index] == 0)
			{
				seed = index;
				seed_counter = i_seed;
				break;
			}
		}
	}
	//����ǩһ����sv�ϲ�
	regions.clear();
	for (map<int, int>::iterator itr = sv_labels_.begin(); itr != sv_labels_.end(); ++itr)
	{
		int label = itr->second;
		int svlabel = itr->first;
		if (label == 0)
			continue;
		regions[label].push_back(sv_clusters_[svlabel]);
		seglabel_to_svlist_[label].push_back(svlabel);
	}
	regions_.erase(regions_.begin(), regions_.end());
	regions_ = regions;
}
//��ȡsegment�������ǩ��ͨ������ÿ��segment�ڲ���sv������label��ͬ��sv���ɱ߽������
void BoundaryVCCS::getRegionNeighbor()
{
	//ѭ��ÿ��region��sv list
	for (map<int, vector<int>>::iterator itr = seglabel_to_svlist_.begin(); itr != seglabel_to_svlist_.end(); ++itr)
	{
		int region_label = itr->first;	
		vector<int> svlist = itr->second;
		for (int i = 0; i < svlist.size(); i++)//ѭ����ǰregion�ڲ�sv
		{
			int svlabel = svlist[i];
			pair<int_multimap, int_multimap> nei_sv_labels = nei_labels_.equal_range(svlabel);//ȷ����ǰsv������
			for (auto iter = nei_sv_labels.first; iter != nei_sv_labels.second; ++iter)//��������
			{
				int cur_sv_label = iter->second;
				int seg_label = sv_labels_[cur_sv_label];//��ȡ����label��Ӧ��region label
				if (seg_label != region_label)//�жϱ�ǩ�Ƿ�һ�£����Ƿ�����
				{
					svp[seg_label]++;
					if (isConvex(svlabel, cur_sv_label))//�ж��Ƿ�Ϊ͹
						convex_svp[seg_label]++;
				}
			}
		}
		//merge͹����
		for (map<int, int>::iterator svp_itr = svp.begin(); svp_itr != svp.end(); ++svp_itr)
		{
			int label = svp_itr->first;
			seg_label_to_neighbor_set_[region_label].insert(label);
			int n1 = svp_itr->second;
			int n2 = convex_svp[label];
			if (double(n2) / n1 >= 0.4)//0.55
			{
				//�ı�ǩ��ɾ������
				//vector<Supervoxel<PointT>::Ptr> cur_svs = regions_[region_label];
				vector<Supervoxel<PointT>::Ptr> nei_svs = regions_[label];
				vector<int> nei_sv_labels = seglabel_to_svlist_[label];
				for (int k = 0; k < nei_svs.size(); k++)
				{
					regions_[region_label].push_back(nei_svs[k]);
					seglabel_to_svlist_[region_label].push_back(nei_sv_labels[k]);
				}
				regions_.erase(label);
				seglabel_to_svlist_.erase(label);
			}
		}
		svp.clear();
		convex_svp.clear();
	}
}
//�ԱȽ�С������Ѱ���ڽ����������кϲ�
void BoundaryVCCS::mergeSmallRegions()
{
	for (map<int, set<int>>::iterator itr = seg_label_to_neighbor_set_.begin(); itr != seg_label_to_neighbor_set_.end(); itr++)
	{
		int label = itr->first;
		set<int> adj_labels = itr->second;
		int max_size = seglabel_to_svlist_[label].size();
		int max_seglabel = label;
		if (seglabel_to_svlist_[label].size() < 3)
		{
			//��������棬���ų�����
			for (set<int>::iterator set_itr = adj_labels.begin(); set_itr != adj_labels.end(); ++set_itr)
			{
				int tmp = seglabel_to_svlist_[*set_itr].size();
				if (tmp > max_size)
				{
					if (tmp > 50)
						continue;
					max_size = tmp;
					max_seglabel = *set_itr;
				}
			}
			if (max_seglabel != label)
			{
				vector<Supervoxel<PointT>::Ptr> nei_svs = regions_[label];
				vector<int> nei_sv_labels = seglabel_to_svlist_[label];
				for (int k = 0; k < nei_svs.size(); k++)
				{
					regions_[max_seglabel].push_back(nei_svs[k]);
					seglabel_to_svlist_[max_seglabel].push_back(nei_sv_labels[k]);
				}
				regions_.erase(label);
				seglabel_to_svlist_.erase(label);
			}
		}
	}
}
//�Ե�ǰ���Ӽ�������������������ص�ǰ�ָ���ڰ���������
int BoundaryVCCS::seedGrowing(int seed, int segment_number, bool flag_c)
{
	queue<int> seeds;
	seeds.push(seed);
	sv_labels_[seed] = segment_number;
	int sv_in_segment = 1;
	pair<int_multimap, int_multimap> seed_nei_labels = nei_labels_.equal_range(seed);
	if (seed_nei_labels.first == seed_nei_labels.second)//û������ֱ�ӷ���
	{
		return sv_in_segment;
	}
	while (!seeds.empty())
	{
		int curr_seed;
		curr_seed = seeds.front();//ȡ������
		seeds.pop();
		//������ǰ��������
		pair<int_multimap, int_multimap> seed_nei_labels = nei_labels_.equal_range(curr_seed);	
		auto iter = seed_nei_labels.first;
		while (iter != seed_nei_labels.second)
		{
			int nei_label = iter->second;
			if (sv_labels_[nei_label] != 0)//�ѱ���ǩ����������ѭ��
			{
				iter++;
				continue;
			}
			bool is_a_seed = false;
			bool in_segment = isValidated(seed, nei_label, is_a_seed, flag_c);
			if (!in_segment)
			{
				iter++;
				continue;
			}
			sv_labels_[nei_label] = segment_number;
			sv_in_segment++;
			if (is_a_seed)
			{
				seeds.push(nei_label);
			}
			iter++;
		}//next neighbor
	}//next seed
	return sv_in_segment;
}
//�ж����������Ƿ���������׼��:����нǣ��д��Ľ�������
bool BoundaryVCCS::isValidated(int seed, int neighbor, bool &is_a_seed, bool flag_c)
{
	is_a_seed = true;
	float cosine_threshold = cos(25.0f / 180.0f * static_cast<float> (M_PI));
	Eigen::Vector3f seed_point = sv_clusters_[seed]->centroid_.getVector3fMap();
	Eigen::Vector3f n_seed = sv_clusters_[seed]->normal_.getNormalVector3fMap();
	Eigen::Vector3f n_nei = sv_clusters_[neighbor]->normal_.getNormalVector3fMap();
	Eigen::Vector3f nei_point = sv_clusters_[neighbor]->centroid_.getVector3fMap();
	//����н�����ֵ,smothness_check and use residual check again pass
	float Dn = abs(n_seed.dot(n_nei));
	if (Dn > cosine_threshold && (seed_point - nei_point).dot(n_seed) < 0.15)
		return true;//< 25��
	//check curture for puching the point into the seed quece
	if (flag_c)
	{
		if (label_curtures_[neighbor] > 0.05)
			is_a_seed = false;
	}
	else
	{
		if (label_wis_[neighbor] > 0.05)
			is_a_seed = false;
	}
	
	return false;
}
//bool BoundaryVCCS::isValidated2(int seed, int neighbor, bool &is_a_seed)
//{
//	is_a_seed = true;
//	Eigen::Vector3f seed_point = sv_clusters_[seed]->centroid_.getVector3fMap();
//	Eigen::Vector3f n_seed = sv_clusters_[seed]->normal_.getNormalVector3fMap();
//	Eigen::Vector3f n_nei = sv_clusters_[neighbor]->normal_.getNormalVector3fMap();
//	Eigen::Vector3f nei_point = sv_clusters_[neighbor]->centroid_.getVector3fMap();
//	Eigen::Vector3f xd;
//	xd = seed_point - nei_point;
//	//smoothness
//	float Dn = pow(n_seed.dot(n_nei), 2);
//	//continuity
//	float Dd = pow((seed_point - nei_point).dot(n_seed), 2) + pow((seed_point - nei_point).dot(n_nei), 2);
//	//convexity
//	float a1 = getAngle3D(xd, n_seed);
//	float a2 = getAngle3D(xd, n_nei);
//	float Dc;
//	if (a1 - a2 > 0)
//	{
//		Dc = (a1 - a2)*(a1 - a2);
//	}
//	else
//	{
//		Dc = (a1 - a2)*(a1 - a2) + (a1 + a2 - PI)*(a1 + a2 - PI);
//	}
//	float k = 0.1;
//	float wi = exp(-Dn / (2 * k*k))*exp(-Dd / (2 * k*k))*exp(-Dc / (2 * k*k));
//	if (wi < 0.001)
//		return true;
//	if (label_wis_[neighbor] > 0.01)
//		is_a_seed = false;
//	return false;
//
//}
//�жϱ߽��ڽ�sv�İ�͹��
bool BoundaryVCCS::isConvex(int s_label, int t_label)
{
	Eigen::Vector3f s_point = sv_clusters_[s_label]->centroid_.getArray3fMap();
	Eigen::Vector3f t_point = sv_clusters_[t_label]->centroid_.getArray3fMap();
	Eigen::Vector3f vec_s_to_t, vec_t_to_s;
	vec_s_to_t = t_point - s_point;
	vec_t_to_s = -vec_s_to_t;	
	Eigen::Vector3f s_normal = sv_clusters_[s_label]->normal_.getNormalVector3fMap();
	Eigen::Vector3f t_normal = sv_clusters_[t_label]->normal_.getNormalVector3fMap();
	Eigen::Vector3f ncross;
	ncross = s_normal.cross(t_normal);
	float normal_angle = getAngle3D(s_normal, t_normal, true);
	//Sanity Criterion
	float intersection_angle = getAngle3D(ncross, vec_t_to_s, true);
	float min_intersect_angle = (intersection_angle < 90.) ? intersection_angle : 180. - intersection_angle;
	float intersect_thresh = 60. * 1. / (1. + exp(-0.25 * (normal_angle - 25.)));
	//Convexity Criterion
	float angle = getAngle3D(vec_t_to_s, s_normal) - getAngle3D(vec_t_to_s, t_normal);
	if (angle <= 0 || normal_angle<10)
		return true;
	/*if (min_intersect_angle>intersect_thresh && (angle <= 0 || normal_angle<10))
		return true;*/
	return false;
}
//��ȡ��ǩ����
PointCloud<PointXYZL>::Ptr BoundaryVCCS::getLabelCloud()
{
	PointCloud<PointXYZL>::Ptr label_cloud(new PointCloud<PointXYZL>);
	for (map<uint32_t, Supervoxel<PointT>::Ptr >::iterator labelsv_itr = sv_clusters_.begin(); labelsv_itr != sv_clusters_.end(); ++labelsv_itr)
	{
		PointCloudT::Ptr voxels = labelsv_itr->second->voxels_;
		int label = labelsv_itr->first;
		PointCloud<PointXYZL> xyzl_copy;
		copyPointCloud(*voxels, xyzl_copy);
		PointCloud<PointXYZL>::iterator xyzl_copy_itr = xyzl_copy.begin();
		for (; xyzl_copy_itr != xyzl_copy.end(); ++xyzl_copy_itr)
			xyzl_copy_itr->label = label;//tag the label for every point in the voxels
		*label_cloud += xyzl_copy;
	}
	return label_cloud;
}
//��ȡÿ��patch/facet��ǩ����
PointCloud<PointXYZL>::Ptr BoundaryVCCS::getPatchCloud()
{
	PointCloud<PointXYZL>::Ptr segment_cloud(new PointCloud<PointXYZL>);
	for (map<int, vector<Supervoxel<PointT>::Ptr>>::iterator itr = regions_.begin(); itr != regions_.end(); ++itr)
	{
		int label = itr->first;
		vector<Supervoxel<PointT>::Ptr> svs = itr->second;	
		PointCloud<PointXYZL> xyzl_copy;
		/*if (svs.size() == 1)
			continue;*/
		for (int i = 0; i < svs.size(); i++)
		{
			PointCloud<PointXYZL> temp;
			copyPointCloud(*svs[i]->voxels_, temp);
			xyzl_copy += temp;
		}
		for (PointCloud<PointXYZL>::iterator xyzl_copy_itr = xyzl_copy.begin(); xyzl_copy_itr != xyzl_copy.end(); ++xyzl_copy_itr)
			xyzl_copy_itr->label = label;//tag the label for every point in the voxels
		*segment_cloud += xyzl_copy;
	}
	return segment_cloud;
}
//��ȡsv��FPFH����
PointCloud<FPFHSignature33>::Ptr BoundaryVCCS::getFPFH(int label)
{
	FPFHEstimation<PointT, Normal, FPFHSignature33> fpfh;
	fpfh.setInputCloud(sv_clusters_[label]->voxels_);
	fpfh.setInputNormals(sv_clusters_[label]->normals_);
	search::KdTree<PointT>::Ptr tree(new search::KdTree<PointT>);
	fpfh.setSearchMethod(tree);
	PointCloud<FPFHSignature33>::Ptr fpfhs(new PointCloud<FPFHSignature33>);
	fpfh.setRadiusSearch(2 * Rvoxel_);
	fpfh.compute(*fpfhs);
	return fpfhs;
}