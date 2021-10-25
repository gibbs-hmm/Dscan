#if !defined (BLOSUM62H)
#define BLOSUM62H

static float blosum62[21][21] = {
/*    C           G           A           S           T       
      N           D           E           Q           K       
      R           H           W           Y           F       
      V           I           L           M           P       */
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,   4.2911,    -1.2502,    -0.2043,    -0.4375,    -0.4333, 
   -1.3299,    -1.7300,    -1.8062,    -1.4509,    -1.5182, 
   -1.6946,    -1.4939,    -1.1521,    -1.2036,    -1.1877, 
   -0.4038,    -0.6138,    -0.6387,    -0.7099,    -1.3976, },
{0,   -1.2502,    2.7816,    0.0798,    -0.1462,    -0.7877, 
   -0.2114,    -0.6568,    -1.0551,    -0.8926,    -0.7640, 
   -1.1521,    -1.0204,    -1.2457,    -1.5199,    -1.5537, 
   -1.5694,    -1.8624,    -1.8135,    -1.3383,    -1.0668, },
{0,   -0.2043,    0.0798,    1.9646,    0.5579,    -0.0227, 
   -0.7654,    -0.8767,    -0.4319,    -0.4020,    -0.3670, 
   -0.7068,    -0.8126,    -1.2634,    -0.8820,    -1.1050, 
   -0.0947,    -0.6609,    -0.7323,    -0.4676,    -0.4071, },
{0,   -0.4375,    -0.1462,    0.5579,    1.9422,    0.6906, 
   0.3005,    -0.1305,    -0.0735,    -0.0506,    -0.1017, 
   -0.3824,    -0.4408,    -1.3759,    -0.8429,    -1.1845, 
   -0.8231,    -1.1741,    -1.2213,    -0.7404,    -0.4045, },
{0,   -0.4333,    -0.7877,    -0.0227,    0.6906,    2.2727, 
   -0.0230,    -0.5254,    -0.4316,    -0.3377,    -0.3348, 
   -0.5612,    -0.8429,    -1.2145,    -0.8030,    -1.0538, 
   -0.0278,    -0.3588,    -0.5987,    -0.3331,    -0.5376, },
{0,   -1.3299,    -0.2114,    -0.7654,    0.3005,    -0.0230, 
   2.8266,    0.6358,    -0.1340,    0.0008,    -0.0895, 
   -0.2199,    0.2892,    -1.8480,    -1.0409,    -1.4970, 
   -1.4382,    -1.6085,    -1.6895,    -1.0754,    -1.0002, },
{0,   -1.7300,    -0.6568,    -0.8767,    -0.1305,    -0.5254, 
   0.6358,    2.8871,    0.7552,    -0.1567,    -0.3509, 
   -0.8029,    -0.5595,    -2.1072,    -1.5325,    -1.7419, 
   -1.5713,    -1.5606,    -1.8028,    -1.5293,    -0.7401, },
{0,   -1.8062,    -1.0551,    -0.4319,    -0.0735,    -0.4316, 
   -0.1340,    0.7552,    2.4514,    0.9273,    0.3877, 
   -0.0577,    -0.0588,    -1.4177,    -1.0102,    -1.5962, 
   -1.2211,    -1.5972,    -1.4232,    -0.9990,    -0.5581, },
{0,   -1.4509,    -0.8926,    -0.4020,    -0.0506,    -0.3377, 
   0.0008,    -0.1567,    0.9273,    2.6426,    0.6363, 
   0.4914,    0.2240,    -0.9732,    -0.7105,    -1.5822, 
   -1.0992,    -1.3848,    -1.0670,    -0.2105,    -0.6410, },
{0,   -1.5182,    -0.7640,    -0.3670,    -0.1017,    -0.3348, 
   -0.0895,    -0.3509,    0.3877,    0.6363,    2.2523, 
   1.0544,    -0.3605,    -1.4782,    -0.9100,    -1.5393, 
   -1.1312,    -1.3351,    -1.2234,    -0.6774,    -0.5068, },
{0,   -1.6946,    -1.1521,    -0.7068,    -0.3824,    -0.5612, 
   -0.2199,    -0.8029,    -0.0577,    0.4914,    1.0544, 
   2.7367,    -0.1249,    -1.3397,    -0.8469,    -1.3932, 
   -1.2513,    -1.4951,    -1.0773,    -0.6836,    -1.0543, },
{0,   -1.4939,    -1.0204,    -0.8126,    -0.4408,    -0.8429, 
   0.2892,    -0.5595,    -0.0588,    0.2240,    -0.3605, 
   -0.1249,    3.7555,    -1.1711,    0.8463,    -0.6171, 
   -1.5587,    -1.6158,    -1.3934,    -0.7756,    -1.0805, },
{0,   -1.1521,    -1.2457,    -1.2634,    -1.3759,    -1.2145, 
   -1.8480,    -2.1072,    -1.4177,    -0.9732,    -1.4782, 
   -1.3397,    -1.1711,    5.2520,    1.0771,    0.4588, 
   -1.4171,    -1.2903,    -0.8159,    -0.7124,    -1.8271, },
{0,   -1.2036,    -1.5199,    -0.8820,    -0.8429,    -0.8030, 
   -1.0409,    -1.5325,    -1.0102,    -0.7105,    -0.9100, 
   -0.8469,    0.8463,    1.0771,    3.2975,    1.4696, 
   -0.6038,    -0.6657,    -0.5310,    -0.4974,    -1.4599, },
{0,   -1.1877,    -1.5537,    -1.1050,    -1.1845,    -1.0538, 
   -1.4970,    -1.7419,    -1.5962,    -1.5822,    -1.5393, 
   -1.3932,    -0.6171,    0.4588,    1.4696,    3.0230, 
   -0.4245,    -0.0804,    0.2074,    0.0063,    -1.7986, },
{0,   -0.4038,    -1.5694,    -0.0947,    -0.8231,    -0.0278, 
   -1.4382,    -1.5713,    -1.2211,    -1.0992,    -1.1312, 
   -1.2513,    -1.5587,    -1.4171,    -0.6038,    -0.4245, 
   1.8845,    1.2735,    0.3942,    0.3436,    -1.1744, },
{0,   -0.6138,    -1.8624,    -0.6609,    -1.1741,    -0.3588, 
   -1.6085,    -1.5606,    -1.5972,    -1.3848,    -1.3351, 
   -1.4951,    -1.6158,    -1.2903,    -0.6657,    -0.0804, 
   1.2735,    1.9993,    0.7608,    0.5634,    -1.3783, },
{0,   -0.6387,    -1.8135,    -0.7323,    -1.2213,    -0.5987, 
   -1.6895,    -1.8028,    -1.4232,    -1.0670,    -1.2234, 
   -1.0773,    -1.3934,    -0.8159,    -0.5310,    0.2074, 
   0.3942,    0.7608,    1.9247,    0.9959,    -1.4300, },
{0,   -0.7099,    -1.3383,    -0.4676,    -0.7404,    -0.3331, 
   -1.0754,    -1.5293,    -0.9990,    -0.2105,    -0.6774, 
   -0.6836,    -0.7756,    -0.7124,    -0.4974,    0.0063, 
   0.3436,    0.5634,    0.9959,    2.6963,    -1.2382, },
{0,   -1.3976,    -1.0668,    -0.4071,    -0.4045,    -0.5376, 
   -1.0002,    -0.7401,    -0.5581,    -0.6410,    -0.5068, 
   -1.0543,    -1.0805,    -1.8271,    -1.4599,    -1.7986, 
   -1.1744,    -1.3783,    -1.4300,    -1.2382,    3.6823, }
};


/* blosum62 = log2(blosum62L)   */

static float blosum62L[21][21] = {
/*    C           G           A           S           T       
      N           D           E           Q           K       
      R           H           W           Y           F       
      V           I           L           M           P       */
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,   19.5772,    0.4204,    0.8680,    0.7384,    0.7406, 
   0.3978,    0.3015,    0.2859,    0.3658,    0.3491, 
   0.3089,    0.3551,    0.4500,    0.4342,    0.4390, 
   0.7559,    0.6535,    0.6423,    0.6114,    0.3796, },
{0,   0.4204,    6.8761,    1.0569,    0.9036,    0.5793, 
   0.8637,    0.6343,    0.4813,    0.5386,    0.5889, 
   0.4500,    0.4930,    0.4217,    0.3487,    0.3406, 
   0.3369,    0.2750,    0.2845,    0.3955,    0.4774, },
{0,   0.8680,    1.0569,    3.9030,    1.4721,    0.9844, 
   0.5883,    0.5446,    0.7413,    0.7568,    0.7754, 
   0.6127,    0.5694,    0.4166,    0.5426,    0.4649, 
   0.9365,    0.6325,    0.6019,    0.7232,    0.7541, },
{0,   0.7384,    0.9036,    1.4721,    3.8429,    1.6140, 
   1.2316,    0.9135,    0.9503,    0.9655,    0.9319, 
   0.7672,    0.7367,    0.3853,    0.5575,    0.4400, 
   0.5652,    0.4432,    0.4289,    0.5986,    0.7555, },
{0,   0.7406,    0.5793,    0.9844,    1.6140,    4.8323, 
   0.9842,    0.6948,    0.7414,    0.7913,    0.7929, 
   0.6777,    0.5575,    0.4309,    0.5732,    0.4817, 
   0.9809,    0.7798,    0.6603,    0.7938,    0.6889, },
{0,   0.3978,    0.8637,    0.5883,    1.2316,    0.9842, 
   7.0940,    1.5538,    0.9113,    1.0006,    0.9398, 
   0.8586,    1.2220,    0.2778,    0.4860,    0.3543, 
   0.3690,    0.3279,    0.3100,    0.4745,    0.4999, },
{0,   0.3015,    0.6343,    0.5446,    0.9135,    0.6948, 
   1.5538,    7.3978,    1.6879,    0.8971,    0.7841, 
   0.5732,    0.6785,    0.2321,    0.3457,    0.2990, 
   0.3365,    0.3390,    0.2866,    0.3464,    0.5987, },
{0,   0.2859,    0.4813,    0.7413,    0.9503,    0.7414, 
   0.9113,    1.6879,    5.4695,    1.9017,    1.3083, 
   0.9608,    0.9601,    0.3743,    0.4965,    0.3307, 
   0.4290,    0.3305,    0.3729,    0.5003,    0.6792, },
{0,   0.3658,    0.5386,    0.7568,    0.9655,    0.7913, 
   1.0006,    0.8971,    1.9017,    6.2446,    1.5543, 
   1.4058,    1.1680,    0.5094,    0.6111,    0.3340, 
   0.4668,    0.3829,    0.4773,    0.8642,    0.6413, },
{0,   0.3491,    0.5889,    0.7754,    0.9319,    0.7929, 
   0.9398,    0.7841,    1.3083,    1.5543,    4.7644, 
   2.0769,    0.7789,    0.3589,    0.5322,    0.3441, 
   0.4565,    0.3964,    0.4283,    0.6253,    0.7038, },
{0,   0.3089,    0.4500,    0.6127,    0.7672,    0.6777, 
   0.8586,    0.5732,    0.9608,    1.4058,    2.0769, 
   6.6654,    0.9171,    0.3951,    0.5560,    0.3807, 
   0.4201,    0.3548,    0.4739,    0.6226,    0.4815, },
{0,   0.3551,    0.4930,    0.5694,    0.7367,    0.5575, 
   1.2220,    0.6785,    0.9601,    1.1680,    0.7789, 
   0.9171,    13.5057,    0.4441,    1.7979,    0.6520, 
   0.3395,    0.3263,    0.3807,    0.5841,    0.4729, },
{0,   0.4500,    0.4217,    0.4166,    0.3853,    0.4309, 
   0.2778,    0.2321,    0.3743,    0.5094,    0.3589, 
   0.3951,    0.4441,    38.1074,    2.1098,    1.3744, 
   0.3745,    0.4089,    0.5681,    0.6103,    0.2818, },
{0,   0.4342,    0.3487,    0.5426,    0.5575,    0.5732, 
   0.4860,    0.3457,    0.4965,    0.6111,    0.5322, 
   0.5560,    1.7979,    2.1098,    9.8321,    2.7695, 
   0.6580,    0.6304,    0.6921,    0.7084,    0.3635, },
{0,   0.4390,    0.3406,    0.4649,    0.4400,    0.4817, 
   0.3543,    0.2990,    0.3307,    0.3340,    0.3441, 
   0.3807,    0.6520,    1.3744,    2.7695,    8.1286, 
   0.7451,    0.9458,    1.1546,    1.0044,    0.2875, },
{0,   0.7559,    0.3369,    0.9365,    0.5652,    0.9809, 
   0.3690,    0.3365,    0.4290,    0.4668,    0.4565, 
   0.4201,    0.3395,    0.3745,    0.6580,    0.7451, 
   3.6922,    2.4175,    1.3142,    1.2689,    0.4431, },
{0,   0.6535,    0.2750,    0.6325,    0.4432,    0.7798, 
   0.3279,    0.3390,    0.3305,    0.3829,    0.3964, 
   0.3548,    0.3263,    0.4089,    0.6304,    0.9458, 
   2.4175,    3.9981,    1.6944,    1.4777,    0.3847, },
{0,   0.6423,    0.2845,    0.6019,    0.4289,    0.6603, 
   0.3100,    0.2866,    0.3729,    0.4773,    0.4283, 
   0.4739,    0.3807,    0.5681,    0.6921,    1.1546, 
   1.3142,    1.6944,    3.7966,    1.9943,    0.3711, },
{0,   0.6114,    0.3955,    0.7232,    0.5986,    0.7938, 
   0.4745,    0.3464,    0.5003,    0.8642,    0.6253, 
   0.6226,    0.5841,    0.6103,    0.7084,    1.0044, 
   1.2689,    1.4777,    1.9943,    6.4814,    0.4239, },
{0,   0.3796,    0.4774,    0.7541,    0.7555,    0.6889, 
   0.4999,    0.5987,    0.6792,    0.6413,    0.7038, 
   0.4815,    0.4729,    0.2818,    0.3635,    0.2875, 
   0.4431,    0.3847,    0.3711,    0.4239,    12.8376, }
};


static float iden[5][5] = 
{
  /* C G A T */
  { 0, 0, 0, 0, 0, },
  { 0, 20, 0, 0, 0, },
  { 0, 0, 20, 0, 0, },
  { 0, 0, 0, 20, 0, },
  { 0, 0, 0, 0, 20, }
};


/* DNA PAM1 mutation probaility matrix - Agarwal and States, J. Comp. Bio. 3 (1), 1-17 1996 */
static double pam1[5][5] =
{
 /* X  C      G      A      T */
  { 0, 0,     0,     0,     0, },
  { 0, 0.99,  0.002, 0.002, 0.006, },
  { 0, 0.002, 0.99,  0.006, 0.002, },
  { 0, 0.002, 0.006, 0.99,  0.002, },
  { 0, 0.006, 0.002, 0.002, 0.99,  }
};


/* DNA PAM1 scoring matrix - Agarwal and States, J. Comp. Bio. 3 (1), 1-17 1996 */
/* this matrix assumes equal background probs - it will be recalculated by */
/* CalcPam1ScoreMatrix for actual background probs in data base */
static double pam1Score[5][5] =
{
 /* X  C      G      A      T */
  { 0, 0,     0,     0,     0, },
  { 0, 1.99,  -6.97, -6.97, -5.38, },
  { 0, -6.97, 1.99,  -5.38, -6.97, },
  { 0, -6.97, -5.38, 1.99,  -6.97, },
  { 0, -5.38, -6.97, -6.97, 1.99,  }
};


/* pam1Score = log2(pam1ScoreL)   */

static double pam1ScoreL[5][5] =
{
 /* X  C      G       A       T */
  { 0, 0,     0,      0,      0, },
  { 0, 3.97,  0.008,  0.008,  0.024, },
  { 0, 0.008, 3.97,   0.024,  0.008, },
  { 0, 0.008, 0.024,  3.97,   0.008, },
  { 0, 0.024, 0.008,  0.008,  3.97,  }
};

#endif