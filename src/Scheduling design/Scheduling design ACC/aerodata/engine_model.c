/*---------------------------------------------------------------------- */
/*                                                                       */
/* F-16 engine model taken from NASA report 1538                         */
/* based on matlab functions by Ying Huo                                 */
/*                                                                       */
/* File "engine_model.c"                                                 */
/* by L. Sonneveldt                                                      */
/* May, 2006                                                             */
/*                                                                       */
/*---------------------------------------------------------------------- */

double tgear (double throt){

    double tgear_value;
    
    if (throt <= 0.77){
    tgear_value = 64.94 * throt;
    }
    else
    {
    tgear_value = 217.38 * throt - 117.38;
    }
    return tgear_value;    
}

double rtau (double dp ){
    
    double tau_r;
    
    if ( dp <= 25.0 ){
        tau_r = 1.0;  
    }
    else if ( dp >= 50.0 ){
        tau_r = 0.1;
    }
    else{
        tau_r = 1.9 - 0.036 * dp;
    }
    return tau_r;
}

double pdot ( double pw, double cpw ){

    double tpw, ta;

    if ( cpw >= 50.0 ){
        if ( pw >= 50.0){
                tpw = cpw;
                ta = 5.0;
            }
        else
        {
            tpw = 60.0;
            ta = rtau ( tpw - pw );
        }
    }
    else{
        if ( pw >= 50.0){
            tpw = 40.0;
            ta = 5.0;
        }
        else
        {
            tpw = cpw;
            ta = rtau ( tpw - pw );
        }
    }
    return ta * ( tpw - pw );
}

double thrst ( double powr, double alti, double rmach ){

    double usalti = alti / 0.3048;

    double Adata[6][6] = {  
    {1060.0,         635.0,          60.0,       -1020.0,       -2700.0,       -3600.0},
    {670.0,         425.0,          25.0,        -710.0,       -1900.0,       -1400.0},
    {880.0,         690.0,         345.0,        -300.0,       -1300.0,        -595.0},
    {1140.0,        1010.0,         755.0,         350.0,        -247.0,        -342.0},
    {1500.0,        1330.0,        1130.0,         910.0,         600.0,        -200.0},
    {1860.0,        1700.0,        1525.0,        1360.0,        1100.0,         700.0}};

    double Bdata[6][6] = {
    {12680.0,       12680.0,       12610.0,       12640.0,       12390.0,       11680.0},
    {9150.0,        9150.0,        9312.0,        9839.0,       10176.0,        9848.0},
    {6200.0,        6313.0,        6610.0,        7090.0,        7750.0,        8050.0},
    {3950.0,        4040.0,        4290.0,        4660.0,        5320.0,        6100.0},
    {2450.0,        2470.0,        2600.0,        2840.0,        3250.0,        3800.0},
    {1400.0,        1400.0,        1560.0,        1660.0,        1930.0,        2310.0}};

    double Cdata[6][6] = {
    {20000.0,       21420.0,       22700.0,       24240.0,       26070.0,       28886.0},
    {15000.0,       15700.0,       16860.0,       18910.0,       21075.0,       23319.0},
    {10800.0,       11225.0,       12250.0,       13760.0,       15975.0,       18300.0},
    {7000.0,        7323.0,        8154.0,        9285.0,       11115.0,       13484.0},
    {4000.0,        4435.0,        5000.0,        5700.0,        6860.0,        8642.0},
    {2500.0,        2600.0,        2835.0,        3215.0,        3950.0,        5057.0}};

int row = 1;
int col = 1;
double h, dh, rm, dm, cdh, s ,t, tmil, tidl, tmax, thrst_value;
int i, m;

h = 0.0001 * usalti;
i = fix( h );
if (i >= 5){
    i = 4;
}
else  if (i <= 0){
     i = 0;
}

dh = h - i;
rm = 5.0 * rmach;

m = fix( rm );
if ( m >= 5 ){
    m = 4;
}
else if (m <= 0){
     m = 0;
}

dm = rm - m;
cdh = 1.0 - dh;

s = Bdata[i+row-1][m+col-1]  * cdh + Bdata[i+1+row-1][m+col-1] * dh;
t = Bdata[i+row-1][m+1+col-1] * cdh + Bdata[i+1+row-1][m+1+col-1] * dh;

tmil = s + ( t - s ) * dm;
if( powr < 50.0 ){
    s = Adata[i+row-1][m+col-1]* cdh + Adata[i+1+row-1][m+col-1] * dh;
    t = Adata[i+row-1][m+1+col-1] * cdh + Adata[i+1+row-1][m+1+col-1] * dh;
    tidl = s + ( t - s ) * dm;
    thrst_value = (tidl + ( tmil - tidl ) * powr *0.02) * 4.4482216 ;
}
else{
    s = Cdata[i+row-1][m+col-1]   * cdh + Cdata[i+1+row-1][m+col-1] * dh;
    t = Cdata[i+row-1][m+1+col-1] * cdh + Cdata[i+1+row-1][m+1+col-1] * dh;
    tmax = s + ( t - s ) * dm;
    thrst_value = (tmil + ( tmax - tmil ) * ( powr - 50.0 ) * 0.02)* 4.4482216;
}
return thrst_value;
}
