/*
allow maximum 512 threads
*/

__kernel void curvature(__global int* l_points, __global int* k_points, 
						__global int* r_point, __global float* couvature, const int size)
{
	const int idx = get_global_id(0);
	if (idx < 512 && 2*idx + 1 < size){
		int i_x = idx;
		while(2*i_x + 1 < size){
			float kl = ((l_points[i_x*2]-k_points[i_x*2])*(l_points[i_x*2]-k_points[i_x*2]) 
						+ (l_points[i_x*2+1]-k_points[i_x*2+1])*(l_points[i_x*2+1]-k_points[i_x*2+1]));
			kl = sqrt(kl);
			float lr = ((l_points[i_x*2]-r_points[i_x*2])*(l_points[i_x*2]-r_points[i_x*2]) 
						+ (l_points[i_x*2+1]-r_points[i_x*2+1])*(l_points[i_x*2+1]-r_points[i_x*2+1]));
			lr = sqrt(lr);
			float kr = ((k_points[i_x*2]-r_points[i_x*2])*(k_points[i_x*2]-r_points[i_x*2]) 
						+ (k_points[i_x*2+1]-r_points[i_x*2+1])*(k_points[i_x*2+1]-r_points[i_x*2+1]));
			kr = sqrt(kr);
			float r = (kl+lr+kr)*(kl+lr-kr)*(kl+kr-lr)*(kr+lr-kl);
			float rayon = kl * lr * kr / sqrt(r);
			int det = (r_points[i_x*2]-k_points[i_x*2])*(l_points[i_x*2+1]-k_points[i_x*2+1])
					-(l_points[i_x*2]-k_points[i_x*2])*(r_points[i_x*2+1]-k_points[i_x*2+1]);
			if (det > 0) {
				couvature[i_x] = 1/rayon;
			}else if (det < 0){
				couvature[i_x] = -1/rayon;
			}else{
				couvature[i_x] = 0.0;
			}

			i_x = i_x + 512;
		}
	}
}