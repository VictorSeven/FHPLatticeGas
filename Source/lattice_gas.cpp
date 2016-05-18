#include<iostream>
#include<cstdlib>
#include<random>
#include<ctime>
#include<fstream>
#include<vector>
#include<string>
#include<chrono>

#define XMAX 5
#define YMAX XMAX*64
#define NMEANDATA 10

using namespace std;

//Inits the variables
int initialize(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> >& nsbit, mt19937 generator,  double prob);

//Collision in the different models
void collisionFHP_I(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator);
void collisionFHP_II(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator);
void collisionFHP_III(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator);

//Get the results for g(d), nu(d), etc for the models
void get_results_FHPI(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);
void get_results_FHPIII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re);

//Particle propagation for any system.
void propagation(vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > cell[7], vector< vector<uint64_t> > nsbit);

//Periodic boundary conditions
int periodic_bc(int k, int maxk);
int periodic_bc(int b, int maxb, int& x);

//Get the k-th bit of n
int bit_at(uint64_t n, int k);

//Write all the data
void writeData(vector< vector<uint64_t> > cell[7], string filename);
void writeGrid();
void writeData(vector<double> data_to_write[7], int N, string filename);
void writeResults(double rho, string filename);

//Measures the most important values of the distribution
void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA]);
void promediate(double mean_data[NMEANDATA], vector<double> ni_to_write[7], vector<double> u_to_write[2], double& rho_eq, double u_eq[2], double Ni_eq[7], int fraction_its, int n_its);

//Gives the equilibrium values for the mean occupation numbers
void equilibrium_ni(double rho, double ci[7][2], double u[2], double ni[7], double nieq[7]);


int main(void)
{
    int i,j;
    double k;
    //Six directions + rest particles
    //  a   b     |   0   1
    //   \ /      |    \ /
    //f --g-- c   |  5--6--2
    //   / \      |    / \
    //  e   d     |   4   3
    vector< vector<uint64_t> > cell[7];
    vector< vector<uint64_t> > result_cell[7];
    vector< vector<uint64_t> > nsbit (XMAX, vector<uint64_t>(YMAX));

    int n_its = 500; //Number of iterations;
    double fraction_its = 0.3; //Fraction of iterations assumed before equilibrium

    //Velocities for the grid above
    double ci[7][2] =
    {
      {-1.0/2.0, sqrt(3.0)/2.0},
      {1.0/2.0, sqrt(3.0)/2.0},
      {1.0, 0.0},
      {1.0/2.0, -sqrt(3.0)/2.0},
      {-1.0/2.0, -sqrt(3.0)/2.0},
      {-1.0, 0.0}
    };

    //Variables to get the numerical data
    double N_particles; //Number of particles

    vector< vector<double> > Ni[7]; //Mean occupation numbers
    double Ni_the[7]; //Theoretical mean occupation numbers in equilibrium

    vector< vector<double> > rho (XMAX*64, vector<double>(YMAX)); //Density
    vector< vector<double> > u[2]; //Velocity
    double mean_data[NMEANDATA]; //Mean of Ni (0-6), rho (7), u (8-9) in all the space
    vector<double> ni_to_write[7]; //Used to write the Ni in every iteration
    vector<double> u_to_write[2]; //Used to write u in every iteration

    double rho_eq; //Values promediated over all the iterations...
    double u_eq[2];
    double Ni_eq[7];

    //Vector init for all the variables
    for (i=0; i < 7; i++)
    {
        cell[i] = vector< vector<uint64_t> >(XMAX, vector<uint64_t>(YMAX));
        result_cell[i] = vector< vector<uint64_t> >(XMAX, vector<uint64_t>(YMAX));
        Ni[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
    }
    for (i=0; i < 2; i++)
    {
        u[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
    }


    //Init of the random generator and distribution
    mt19937 gen (chrono::high_resolution_clock::now().time_since_epoch().count());

    k = 0.2;

    //for (k=0.01; k < 0.5; k += 0.02)
    //{
        //Init the grid and get all the particles
        N_particles = initialize(cell, result_cell, nsbit, gen, k);

        //Init values for promediating for every iteration
        rho_eq = 0.0;
        u_eq[0] = 0.0;
        u_eq[1] = 0.0;
        for (i=0; i <7;i++)
        {
            Ni_eq[i] = 0.0;
        }

        //Do iterations
        for (i=0; i <= n_its; i++)
        {
            collisionFHP_I(cell, result_cell, nsbit, gen); //Collision
            propagation(result_cell,cell,nsbit); //Propagation

            measure(cell, ci, Ni, rho, u, mean_data); //Measure
            promediate(mean_data, ni_to_write, u_to_write, rho_eq, u_eq, Ni_eq, fraction_its, n_its); //Promediate values to obtain _eq ones

            //To do animations
            //writeData(cell,"data"+to_string(i)+".txt");

        }

        //Write ni in every iteration
        //writeData(ni_to_write, 7, "ni");

        //Write u in every iteration
        writeData(u_to_write, 2, "speed");


        //Finish the mean values
        for (j=0; j < 7; j++)
        {
            Ni_eq[j] /= (1.0-fraction_its) * n_its;
        }
        rho_eq /= (1.0-fraction_its) * n_its;
        u_eq[0] /= (1.0-fraction_its) * n_its;
        u_eq[1] /= (1.0-fraction_its) * n_its;


        cout << rho_eq << endl;

        //Calculate equilibrium values for the numerical data obtained
        equilibrium_ni(rho_eq, ci, u_eq, Ni_eq, Ni_the);

        //Append the results to a file.
        //writeResults(rho_eq, "results3.txt");
    //}

    return 0;
}

//Initializces the cells with random speed
int initialize(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> >& nsbit, mt19937 generator, double prob)
{
    int i,j,k,b; //Counters
    int n_particles;
    uint64_t aux_bit; //Auxiliar variables
    uint64_t bit_to_add;
    uniform_real_distribution<double> rbit(0,1); //Used to add random bits

    n_particles = 0; //Counter of particles
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            for (k=0; k < 7; k++)
            {
                uint64_t aux_bit = 0; //Init
                //Get a random number with a p% of ones
                for (b=0; b < 64; b++)
                {
                    bit_to_add = rbit(generator) <= prob ? uint64_t(1) : uint64_t(0);
                    aux_bit = bit_to_add ^ (aux_bit << 1); //Add the one or zero
                }
                cell[k][i][j] = aux_bit; //Add to the cell
                n_particles += __builtin_popcount(cell[k][i][j]); //Count number of ones
                result_cell[k][i][j] = 0;
            }

            nsbit[i][j] = ~uint64_t(0); //Non-solid = Non-wall -> 1 in the fluid
        }
    }
    return n_particles;
}

//Does the FHP-I model collision: double and triple particles.
//Uses the data in the cell, computes collisions and stores it in result_cell.
void collisionFHP_I(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator)
{
    int x,y;

    uint64_t rnd, no_rnd, nsb; //Random, negate of random, and non-solid bit

    uint64_t a,b,c,d,e,f; //Single cells:
    uint64_t tb_ad, tb_be, tb_cf; //Two-body possible collisions
    uint64_t triple; //Three body collision
    uint64_t cha, chb, chc, chd, che, chf; //Change in each cell

    uniform_int_distribution<uint64_t> dist(uint64_t(0), ~uint64_t(0)); //from 000...0000 to 111...1111

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {

            //Shorter names
            a = cell[0][x][y];
            b = cell[1][x][y];
            c = cell[2][x][y];
            d = cell[3][x][y];
            e = cell[4][x][y];
            f = cell[5][x][y];

            //If there're particles in i and i+3 but not in any other cell, then
            //we have a two body collision
            tb_ad = (a&d)&(~(b|c|e|f));
            tb_be = (b&e)&(~(a|c|d|f));
            tb_cf = (c&f)&(~(a|b|d|e));


            //If we don't have any contiguous occupied or non-occupied sites, then
            //we have a three body collision
            triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

            //Get the random bit and wall bit
            rnd = dist(generator);
            no_rnd = ~rnd; //Save time!
            nsb = nsbit[x][y];


            //Change the configuration using the collisions.
            cha = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)&nsb);
            chd = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)&nsb);
            chb = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)&nsb);
            che = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)&nsb);
            chc = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)&nsb);
            chf = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)&nsb);

            //Update the configuration with the change
            result_cell[0][x][y] = a ^ cha;
            result_cell[1][x][y] = b ^ chb;
            result_cell[2][x][y] = c ^ chc;
            result_cell[3][x][y] = d ^ chd;
            result_cell[4][x][y] = e ^ che;
            result_cell[5][x][y] = f ^ chf;
        }
    }

    return;
}

void collisionFHP_II(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator)
{
    int x,y;

    uint64_t rnd, no_rnd, nsb; //Random, negate of random, and non-solid bit

    uint64_t a,b,c,d,e,f,r; //Single cells:
    uint64_t tb_ad, tb_be, tb_cf, tb_ra, tb_rb, tb_rc, tb_rd, tb_re, tb_rf; //Two-body possible collisions
    uint64_t ra, rb, rc, rd, re, rf; //Collision of particles giving a rest one
    uint64_t triple; //Three body collision
    uint64_t cha, chb, chc, chd, che, chf, chr; //Change in each cell

    uniform_int_distribution<uint64_t> dist(uint64_t(0), ~uint64_t(0)); //from 000...0000 to 111...1111

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {

            //Shorter names
            a = cell[0][x][y];
            b = cell[1][x][y];
            c = cell[2][x][y];
            d = cell[3][x][y];
            e = cell[4][x][y];
            f = cell[5][x][y];
            r = cell[6][x][y];

            //If there're particles in i and i+3 but not in any other cell, then
            //we have a two body collision
            tb_ad = (a&d)&(~(b|c|e|f));
            tb_be = (b&e)&(~(a|c|d|f));
            tb_cf = (c&f)&(~(a|b|d|e));
            //With rest:
            tb_ra = (r&a&~(b|c|d|e|f));
            tb_rb = (r&b&~(a|c|d|e|f));
            tb_rc = (r&c&~(a|b|d|e|f));
            tb_rd = (r&d&~(a|b|c|e|f));
            tb_re = (r&e&~(a|b|c|d|f));
            tb_rf = (r&f&~(a|b|c|d|e));


            //If we don't have any contiguous occupied or non-occupied sites, then
            //we have a three body collision
            triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

            //Collision of a two particles at i and i+2 which leads to a rest and a moving one:
            ra = (f&b&~(r|a|c|d|e));
            rb = (a&c&~(r|b|d|e|f));
            rc = (b&d&~(r|a|c|e|f));
            rd = (c&e&~(r|a|b|d|f));
            re = (d&f&~(r|a|b|c|e));
            rf = (e&a&~(r|b|c|d|f));


            //Get the random bit and wall bit
            rnd = dist(generator);
            no_rnd = ~rnd; //Save time!
            nsb = nsbit[x][y];


            //Change the configuration using the collisions.
            cha = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_ra|tb_rb|tb_rf|ra|rb|rf&nsb);
            chd = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_rd|tb_rc|tb_re|rd|rc|re&nsb);
            chb = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_rb|tb_ra|tb_rc|rb|ra|rc&nsb);
            che = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_re|tb_rd|tb_rf|re|rd|rf&nsb);
            chc = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rc|tb_rb|tb_rd|rc|rb|rd&nsb);
            chf = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rf|tb_ra|tb_re|rf|ra|re&nsb);
            chr = (tb_ra|tb_rb|tb_rc|tb_rd|tb_re|tb_rf|ra|rb|rc|rd|re|rf);

            //Update the configuration with the change
            result_cell[0][x][y] = a ^ cha;
            result_cell[1][x][y] = b ^ chb;
            result_cell[2][x][y] = c ^ chc;
            result_cell[3][x][y] = d ^ chd;
            result_cell[4][x][y] = e ^ che;
            result_cell[5][x][y] = f ^ chf;
            result_cell[6][x][y] = r ^ chr;
        }
    }

    return;
}

void collisionFHP_III(vector< vector<uint64_t> > cell[7],  vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > nsbit, mt19937 generator)
{
    int x,y;

    uint64_t rnd, no_rnd, nsb; //Random, negate of random, and non-solid bit

    uint64_t a,b,c,d,e,f,r; //Single cells:
    uint64_t ad, be, cf; //Auxiliar cells;
    uint64_t tb_ad, tb_be, tb_cf, tb_ra, tb_rb, tb_rc, tb_rd, tb_re, tb_rf; //Two-body possible collisions
    uint64_t fb_adbe, fb_adcf, fb_becf; //Possible four-body collisions
    uint64_t ra, rb, rc, rd, re, rf; //Collision of particles giving a rest one
    uint64_t spect_ad_b, spect_ad_c, spect_ad_e, spect_ad_f;
    uint64_t spect_be_a, spect_be_c, spect_be_d, spect_be_f;
    uint64_t spect_cf_a, spect_cf_b, spect_cf_d, spect_cf_e; //Two body head-on collision with spectator
    uint64_t spectator_change_ad, spectator_change_be, spectator_change_cf; //Changes due to spectators
    uint64_t fb_change_ad, fb_change_be, fb_change_cf; //Changes due to 4-body collisions
    uint64_t triple; //Three body collision
    uint64_t cha, chb, chc, chd, che, chf, chr; //Change in each cell

    uniform_int_distribution<uint64_t> dist(uint64_t(0), ~uint64_t(0)); //from 000...0000 to 111...1111

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {
             //Shorter names
            a = cell[0][x][y];
            b = cell[1][x][y];
            c = cell[2][x][y];
            d = cell[3][x][y];
            e = cell[4][x][y];
            f = cell[5][x][y];
            r = cell[6][x][y];

            ad = a&d;
            be = b&e;
            cf = c&f;

            //Get the random bit and wall bit
            rnd = dist(generator);
            no_rnd = ~rnd; //Save time!
            nsb = nsbit[x][y];

            //If there're particles in i and i+3 but not in any other cell, then
            //we have a two body collision
            tb_ad = ad&(~(b|c|e|f));
            tb_be = be&(~(a|c|d|f));
            tb_cf = cf&(~(a|b|d|e));

            //With rest:
            tb_ra = (r&a&~(b|c|d|e|f));
            tb_rb = (r&b&~(a|c|d|e|f));
            tb_rc = (r&c&~(a|b|d|e|f));
            tb_rd = (r&d&~(a|b|c|e|f));
            tb_re = (r&e&~(a|b|c|d|f));
            tb_rf = (r&f&~(a|b|c|d|e));

            //Handle the four body collisions
            fb_adbe = ad&be&(~(c|f));
            fb_adcf = ad&cf&(~(b|e));
            fb_becf = be&cf&(~(a|d));

            //Changes to do due to this collisions. When the ad is present, it changes with p=0.5. In becf it always appears:
            fb_change_ad = (rnd&fb_adbe)|(no_rnd&fb_adcf)|fb_becf;
            fb_change_be = (rnd&fb_becf)|(no_rnd&fb_adbe)|fb_adcf;
            fb_change_cf = (rnd&fb_adcf)|(no_rnd&fb_becf)|fb_adbe;

            //Three body collisions with spectator. The array-form is for legibility. Each two body collision can have 4-subcases.
            spect_ad_b = ad&b&(~(c|e|f));
            spect_ad_c = ad&c&(~(b|e|f));
            spect_ad_e = ad&e&(~(b|c|f));
            spect_ad_f = ad&f&(~(b|c|e));
            spect_be_a = be&a&(~(c|d|f));
            spect_be_c = be&c&(~(a|d|f));
            spect_be_d = be&d&(~(a|c|f));
            spect_be_f = be&f&(~(a|c|d));
            spect_cf_a = cf&a&(~(b|d|e));
            spect_cf_b = cf&b&(~(a|d|e));
            spect_cf_d = cf&d&(~(a|b|e));
            spect_cf_e = cf&e&(~(a|b|d));

            //Changes to do due to this collisions. The "a" will change if we have a head-on collision in ad, or if we have a head on collision
            //in be or cf which leads to a rotation to ad
            spectator_change_ad = (spect_ad_b|spect_ad_e|spect_ad_c|spect_ad_f)|(spect_be_c|spect_be_f)|(spect_cf_b|spect_cf_e);
            spectator_change_be = (spect_be_a|spect_be_c|spect_be_d|spect_be_f)|(spect_ad_c|spect_ad_f)|(spect_cf_a|spect_cf_d);
            spectator_change_cf = (spect_cf_a|spect_cf_b|spect_cf_d|spect_cf_e)|(spect_ad_b|spect_ad_e)|(spect_be_a|spect_be_d);

            //If we don't have any contiguous occupied or non-occupied sites, then
            //we have a three body collision
            triple = (a^b)&(b^c)&(c^d)&(d^e)&(e^f);

            //Collision of a two particles at i and i+2 which leads to a rest and a moving one:
            ra = (f&b&~(r|a|c|d|e));
            rb = (a&c&~(r|b|d|e|f));
            rc = (b&d&~(r|a|c|e|f));
            rd = (c&e&~(r|a|b|d|f));
            re = (d&f&~(r|a|b|c|e));
            rf = (e&a&~(r|b|c|d|f));


            //Change the configuration using the collisions.
            cha = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_ra|tb_rb|tb_rf|ra|rb|rf|spectator_change_ad|fb_change_ad&nsb);
            chd = (tb_ad|triple|(rnd&tb_be)|(no_rnd&tb_cf)|tb_rd|tb_rc|tb_re|rd|rc|re|spectator_change_ad|fb_change_ad&nsb);
            chb = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_rb|tb_ra|tb_rc|rb|ra|rc|spectator_change_be|fb_change_be&nsb);
            che = (tb_be|triple|(rnd&tb_cf)|(no_rnd&tb_ad)|tb_re|tb_rd|tb_rf|re|rd|rf|spectator_change_be|fb_change_be&nsb);
            chc = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rc|tb_rb|tb_rd|rc|rb|rd|spectator_change_cf|fb_change_cf&nsb);
            chf = (tb_cf|triple|(rnd&tb_ad)|(no_rnd&tb_be)|tb_rf|tb_ra|tb_re|rf|ra|re|spectator_change_cf|fb_change_cf&nsb);
            chr = (tb_ra|tb_rb|tb_rc|tb_rd|tb_re|tb_rf|ra|rb|rc|rd|re|rf);

            //Update the configuration with the change
            result_cell[0][x][y] = a ^ cha;
            result_cell[1][x][y] = b ^ chb;
            result_cell[2][x][y] = c ^ chc;
            result_cell[3][x][y] = d ^ chd;
            result_cell[4][x][y] = e ^ che;
            result_cell[5][x][y] = f ^ chf;
            result_cell[6][x][y] = r ^ chr;
        }
    }
}




//Propagation of particles. Result is exported again to cell from the result_cell
//computed in collisions.
void propagation(vector< vector<uint64_t> > result_cell[7], vector< vector<uint64_t> > cell[7], vector< vector<uint64_t> > nsbit)
{
    int x,y, xnext, xprev, ynext, yprev;

    for (x=0; x < XMAX; x++)
    {
        for (y=0; y < YMAX; y++)
        {
            //Use periodic boundary conditions
            xnext = periodic_bc(x+1, XMAX);
            ynext = periodic_bc(y+1, YMAX);
            xprev = periodic_bc(x-1, XMAX);
            yprev = periodic_bc(y-1, YMAX);

            //Propagator:
            cell[0][x][ynext] = (result_cell[0][x][y]<<1)^(result_cell[0][xnext][y]>>63);
            cell[1][x][ynext] = result_cell[1][x][y];
            cell[2][xnext][y] = (result_cell[2][x][y]<<63)^(result_cell[2][xnext][y]>>1);
            cell[3][x][yprev] = (result_cell[3][x][y]>>1)^(result_cell[3][xprev][y]<<63);
            cell[4][x][yprev] = result_cell[4][x][y];
            cell[5][xprev][y] = (result_cell[5][x][y]>>63)^(result_cell[5][xprev][y]<<1);
            cell[6][x][y] = result_cell[6][x][y]; //Rest particles don't move

        }
    }
}

//Computes periodic boundary conditions.
int periodic_bc(int k, int maxk)
{
    if (k < 0)
    {
        return maxk-1;
    }
    else if(k == maxk)
    {
        return 0;
    }
    else
    {
        return k;
    }
}

//Periodic conditions for bit-by-bit checking
int periodic_bc(int b, int maxb, int& x)
{
    if (b < 0)
    {
        x = periodic_bc(x-1, XMAX); //Goes backward on x
        return maxb-1; //Get last bit
    }
    else if(b == maxb)
    {
        x = periodic_bc(x+1, XMAX); //Goes forward on x
        return 0; //Get first bit
    }
    else
    {
        return b;
    }
}

void writeGrid()
{
    ofstream output;

    double L = 1.0; //side of the grid
    int i,j,k; //counters
    double xpos; //Position to write.

    output.open("grid.txt");
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            for (k=0; k < 64; k++)
            {
                xpos = k + 64 * i; //Get every position inside the bit array
                output << 0.5 * j * L + xpos * L << " " << j * L * 0.8660254038 << endl; //Write the output
            }
        }
    }
    output.close();

    return;
}

void writeData(vector< vector<uint64_t> > cell[7], string filename)
{
    ofstream output;
    uint64_t particle; //Is there a particle in?
    double L = 1.0; //side of the grid
    int i,j,k; //counters
    int b; //Every bit of particle
    double xpos; //Position to write.

    output.open(filename);
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {
            //Particle = 1 if there's at least one particle
            particle = cell[0][i][j]|cell[1][i][j]|cell[2][i][j]|cell[3][i][j]|cell[4][i][j]|cell[5][i][j]|cell[6][i][j];

            for (k=0; k < 64; k++)
            {
                b = bit_at(particle, k);
                if (b == 1)
                {
                    xpos = k + 64 * i; //Get every position inside the bit array
                    output << 0.5 * j * L + xpos * L << " " << j * L * 0.8660254038 << endl; //Write the output
                }
            }
        }
    }
    output.close();


    return;
}

//Write the first N data stored in the vector. With N=7, Ni can be written, and
//with N=2, also u
void writeData(vector<double> data_to_write[7], int N, string filename)
{
    ofstream output;
    int i,j;


    for (i=0; i < N; i++)
    {
        output.open(filename + to_string(i) + ".txt");
        for (j=0; j < data_to_write[i].size(); j++)
        {
            output << j << " " << data_to_write[i][j] << endl;
        }
        output.close();
    }


    return;
}

void writeResults(double rho, string filename)
{
    double d, cs, gd, nud, etad, re;
    ofstream file;

    get_results_FHPIII(rho, d, cs, gd, nud, etad, re);

    file.open(filename, ios_base::app);
    file << d << " " << cs << " " << gd << " " << nud << " " << etad << " " << re << endl;
    file.close();

    return;
}

//Returns the k-th bit of n. What it does is to put the k-th as first,
//then check the k-th & 1 (because 1<<63 is 1000...000) and we can obtain then
//0000...000 or 1000...0000 so we do again >>63 to get 1 or 0.
int bit_at(uint64_t n, int k)
{
    return (((n<<k) & ((uint64_t(1)<<63))) >>63);
}

//Measures the mean occupation numbers, density and velocity
void measure(vector< vector<uint64_t> > cell[7], double ci[7][2], vector< vector<double> > Ni[7], vector< vector<double> >&  rho, vector< vector<double> >  u[2], double mean_data[NMEANDATA])
{
    int i,j,k,b;
    int bnext, bprev, ynext, yprev;
    int neiga, neigb, neigc, neigd, neige, neigf; //Following the directions of the FHP grid
    double sum;
    vector< vector<double> > flux[2];

    for (i=0; i < 2; i++)
    {
        flux[i] = vector< vector<double> >(XMAX*64, vector<double>(YMAX));
    }

    //Init sums for the mean data
    for (i=0; i < NMEANDATA; i++)
    {
            mean_data[i] = 0.0;
    }
    //Make measures for every (x,y) in space
    for (i=0; i < XMAX; i++)
    {
        for (j=0; j < YMAX; j++)
        {

            for (b=0; b < 64; b++)
            {
                rho[b + i*64][j] = 0.0; //Initalize the rho in this  (x,y)
                flux[0][b + i*64][j] = 0.0;
                flux[1][b + i*64][j] = 0.0; //The same for flux vector

                bnext = periodic_bc(b+1, 64, i);
                bprev = periodic_bc(b-1, 64, i);
                ynext = periodic_bc(j+1, YMAX);
                yprev = periodic_bc(j-1, YMAX);

                //Calculate Ni for every direction:
                for (k=0; k < 7; k++)
                {
                    neiga = bit_at(cell[k][i][ynext],bprev);
                    neigb = bit_at(cell[k][i][ynext],b);
                    neigc = bit_at(cell[k][i][j],bnext);
                    neigd = bit_at(cell[k][i][yprev],bnext);
                    neige = bit_at(cell[k][i][yprev],b);
                    neigf = bit_at(cell[k][i][j],bprev);

                    Ni[k][b + i*64][j] = (1.0*(neiga+neigb+neigc+neigd+neige+neigf))/6.0; //Obtain the Ni at x,y
                    rho[b + i*64][j]  += Ni[k][b + i*64][j]; //Add the result to rho

                    mean_data[k] += Ni[k][b + i*64][j]; //Computes sum of all Ni for the entire space

                    flux[0][b + i*64][j] += ci[k][0] * Ni[k][b + i*64][j];
                    flux[1][b + i*64][j] += ci[k][1] * Ni[k][b + i*64][j]; //The same for the vector flux

                }
                //Get the velocities from the data calculated, if density not equals 0.0
                u[0][b + i*64][j] = rho[b + i*64][j] != 0 ? flux[0][b + i*64][j] / rho[b + i*64][j] : 0.0;
                u[1][b + i*64][j] = rho[b + i*64][j] != 0 ? flux[1][b + i*64][j] / rho[b + i*64][j] : 0.0;



                //Sum for the entire space for rho and u
                mean_data[7] += rho[b + i*64][j];
                mean_data[8] += u[0][b + i*64][j];
                mean_data[9] += u[1][b + i*64][j];


            }


        }
    }

    //Finish the mean dividing by the number of cells
    for (i=0; i < NMEANDATA; i++)
    {
            mean_data[i] /= (XMAX*64)*YMAX;
    }

    return;
}

//Promediate over measured values in the main loop to get stable data
void promediate(double mean_data[NMEANDATA], vector<double> ni_to_write[7], vector<double> u_to_write[2], double& rho_eq, double u_eq[2], double Ni_eq[7], int fraction_its, int n_its)
{
    int i,j;

    //Promediate in the equilibrium
    if (i >= fraction_its * n_its)
    {
        for (j=0; j < 7; j++)
        {
            Ni_eq[j] += mean_data[j];
        }
        rho_eq += mean_data[7];
        u_eq[0] += mean_data[8];
        u_eq[1] += mean_data[9];
    }

    //Division by the number of particles will be done after the main loop, to save time

    //Add the ni and u measured values to do a graph
    for (j=0; j < 7; j++)
    {
        ni_to_write[j].push_back(mean_data[j]);
    }
    u_to_write[0].push_back(mean_data[8]);
    u_to_write[1].push_back(mean_data[9]);

    return;
}

    //Writes the values on the variables given by reference.
    void get_results_FHPI(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
    {
        d = rho / 6.0;
        cs = 1.0/ sqrt(2.0);
        gd = ((1.0 - 2.0 * d) / (1.0 - d)) / 2.0;
        nud = (1.0 / (d *  pow(1.0 - d, 3.0))) / 12.0 - 1.0/8.0;
        etad = 0.0;
        re = cs * gd / nud;

        return;
    }

    void get_results_FHPII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
    {
        d = rho / 7.0;
        cs = sqrt(3.0 / 7.0);
        gd = 7.0*(1.0-2.0*d)/(12.0 * (1.0 - d));
        nud = 1.0 / (28.0 * d * pow(1.0 - d, 3.0) * (1.0 - 4.0 * d / 7.0)) - 1.0/8.0;
        etad = 1.0 / (98.0 * d * pow(1.0 - d, 4.0)) - 1.0/28.0;
        re = cs * gd / nud;

        return;
    }

    void get_results_FHPIII(double rho, double& d, double& cs, double& gd, double& nud, double& etad, double& re)
    {
        d = rho / 7.0;
        cs = sqrt(3.0 / 7.0);
        gd = 7.0*(1.0-2.0*d)/(12.0 * (1.0 - d));
        nud = 1.0 / (28.0 * d * (1.0 - d) * (1.0 - 8.0 * d * (1.0 - d)/ 7.0)) - 1.0/8.0;
        etad = 1.0 / (98.0 * d * (1.0 - d) * (1.0 - 2.0 * d * (1.0 - d))) - 1.0/28.0;
        re = cs * gd / nud;

        return;
    }


void equilibrium_ni(double rho, double ci[7][2], double u[2], double ni[7], double nieq[7])
{
    int i;
    double Grho = (6.0-2.0*rho)/(3.0*(6-rho));
    double cdotu; //Scalar product ci * u
    for (i=0; i < 7; i++)
    {
        cdotu = ci[i][0]*u[0]+ci[i][1]*u[1];
        nieq[i] = rho/6.0 + rho*(cdotu)/3.0 + rho*Grho*((cdotu*cdotu) - (u[0]*u[0]+u[1]*u[1])/2.0);
    }
    return;
