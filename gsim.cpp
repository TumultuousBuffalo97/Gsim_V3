#include <iostream>
#include <random>
#include <chrono>
#include <ctime>
#include <thread>
#ifdef SP_MATH
#define double float
#endif
typedef struct part {
	double xp, yp, xv, yv;
} part;
void part_calcForce(part* p1, part* p2, double dt) {
	double dx, dy, f;
	dx = p2->xp - p1->xp;
	dy = p2->yp - p1->yp;
	f = (dx * dx + dy * dy) * dt;
	p1->xv += dx / f;
	p1->yv += dy / f;
	p2->xv -= dx / f;
	p2->yv -= dy / f;
}
void part_calcForce_mt(part* parts, int nparts, int pstart, int pend, double dt) {
	int p1, p2;
	for (p1 = pstart; p1 < pend; ++p1) {
		for (p2 = p1 + 1; p2 < nparts; ++p2) {
			part_calcForce(&parts[p1], &parts[p2], dt);
		}
	}
}
int main(int argc, char* argv[]) {
	printf("Initializing...\n");
	FILE* fp = fopen("./out.gsim", "wb");//open output file
	//read arguments from command line
	int nframes = atoi(argv[1]);
	int xr = atoi(argv[2]);;
	int yr = atoi(argv[3]);;
	int nparts = atoi(argv[4]);;
	double dt = atof(argv[5]);
	int nthreads = atoi(argv[6]);
	dt *= 60.0;
	//write sim info into file
	fwrite(&nframes, sizeof(int), 1, fp);
	fwrite(&xr, sizeof(int), 1, fp);
	fwrite(&yr, sizeof(int), 1, fp);
	fwrite(&nparts, sizeof(int), 1, fp);

	part* parts = (part*)calloc(nparts, sizeof(part));//allocate memory for particles
	{//initialize parts
		//make a random number generator
		std::default_random_engine random;
		std::uniform_real_distribution<double> distribution(0, 1);
		
		int p;
		printf("Assigning particle positions...\n");
		for (p = 0; p < nparts; ++p) {//set particles randomly
			parts[p].xp = distribution(random) * ((double)xr - 1.);
			parts[p].yp = distribution(random) * ((double)yr - 1.);
			parts[p].xv = 0;
			parts[p].yv = 0;
		}
	}
	printf("Scheduling threads...\n");
	std::thread* threads = (std::thread*)calloc(nthreads, sizeof(std::thread));//threads
	int* tRanges = (int*)calloc(sizeof(int), nthreads + 1);//particle start and end for each thread
	{//find particle ranges for threads
		//find what particle range would give an evenly divided number of interations between each thread
		unsigned long long nInters = 0.5 * nparts * nparts - 0.5 * nparts;
		unsigned long long nIntersPerThread = nInters / nthreads;
		unsigned long long count = 0;
		int ntd = 1;
		int p1, p2;
		if (nthreads>1) {
			for (p1 = 0; p1 < nparts; ++p1) {
				for (p2 = p1 + 1; p2 < nparts; ++p2) {
					count++;
					if (count > ntd * nIntersPerThread) {
						tRanges[ntd] = p1;
						ntd++;
					}
				}
			}
			tRanges[ntd] = nparts;
		}
		else {
			tRanges[0] = 0;
			tRanges[1] = nparts;
		}
	}
	printf("Simulating...\n");
	{
		//clock stuff
		auto start = std::chrono::steady_clock::now();
		auto end = std::chrono::steady_clock::now();
		std::chrono::milliseconds elapsed;

		int f;//current frame
		int p;//working particle
		int t;//thread index
		double nx, ny;//next x and y for detecting edge collisions
		float xpf, ypf;//used for converting positions from double to single precision
		for (f = 0; f < nframes; ++f) {
			start = std::chrono::steady_clock::now();//get time at beginning of frame

			for (t = 0; t < nthreads - 1; ++t) {//run each thread
				threads[t] = std::thread(part_calcForce_mt, parts, nparts, tRanges[t], tRanges[t + 1], dt);
			}
			threads[t] = std::thread(part_calcForce_mt, parts, nparts, tRanges[t], tRanges[t + 1], dt);
			for (t = 0; t < nthreads; ++t) {//wait for each thread to finish
				threads[t].join();
			}

			for (p = 0; p < nparts; ++p) {//loop through particles
				//check if the particles are going to be out of bounds
				nx = parts[p].xp + parts[p].xv;
				ny = parts[p].yp + parts[p].yv;
				//reflect them off the bounds
				if (nx >= xr | nx <= 0) {
					parts[p].xv *= -.5;
				}
				if (ny >= yr | ny <= 0) {
					parts[p].yv *= -.5;
				}
				//convert to float from double and write to file
				xpf = parts[p].xp;
				ypf = parts[p].yp;
				fwrite(&xpf, sizeof(float), 1, fp);
				fwrite(&ypf, sizeof(float), 1, fp);
				
				//advance particles forward
				parts[p].xp += parts[p].xv;
				parts[p].yp += parts[p].yv;
			}
			
			end = std::chrono::steady_clock::now();
			elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			printf("frame:%d, frametime:%llums, framerate: %f\r", f, elapsed.count(), 1000.0 / elapsed.count());
		}
	}
	return 0;
}