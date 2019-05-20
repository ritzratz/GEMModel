/* ===============================================================*/
/* GEM efficiency calculations v3.0                               */
/*                                                                */
/* Written by Viktor Ratza, ratza@hiskp.uni-bonn.de, May 2019     */
/* ===============================================================*/

class GEM_Efficiencies_v3
{
    public:

      GEM_Efficiencies_v3(int _N);

      void   SetGeometry(double _d, double _p, double _L, double _g1, double _g2);
      void   SetN(int _N);

      double GetParameter_C1();
      double GetParameter_C2();
      double GetParameter_C3(); 
      double GetParameter_C1Bar();
      double GetParameter_C2Bar();
      double GetParameter_C3Bar();
      double GetParameter_C4();
      double GetParameter_C5();
      double GetParameter_C6();
      double GetParameter_C7(double eta1, double eta2); 
      double GetParameter_C8(double eta1, double eta2); 
      double GetParameter_C9(double eta1, double eta2);     
      double GetParameter_C7Bar_fromx(double xa, double xb);
      double GetParameter_C8Bar_fromx(double xa, double xb);
      double GetParameter_C9Bar_fromx(double xa, double xb);
      double GetParameter_C7Bar(double eta1, double eta2);  
      double GetParameter_C8Bar(double eta1, double eta2);
      double GetParameter_C9Bar(double eta1, double eta2);
       
      double Get_xb_top(double eta1, double eta2); 
      double Get_eta1k1(double eta2);
      double Get_eta1k2(double eta2);
      double Get_xb_bot(double eta1, double eta2); 
      double Get_eta2k1(double eta1);
      double Get_eta2k2(double eta1);
      
      double Get_CollTop(double eta1, double eta2); //Scan eta1
      double Get_ExtrBot(double eta1, double eta2); //Scan eta2
      double Get_CollBot(double eta1, double eta2); //Scan eta2
      double Get_ExtrTop(double eta1, double eta2); //Scan eta1
      double Get_Coll_Driftions_g(double eta1, double eta2, double x1, double x2);
      
      //General at cathode 
      double fg_cath_lambda(double xa, double xb);
      double fg_cath_mu1(double xa, double xb);
      double fg_cath_mu2(double xa, double xb);
      double Fg_cath_lambda(int n, double xa, double xb);
      double Fg_cath_mu1(int n, double xa, double xb);
      double Fg_cath_mu2(int n, double xa, double xb);
      double lambda_cath_g(double xa, double xb);
      double mu1_cath_g(double xa, double xb);
      double mu2_cath_g(double xa, double xb);
      
      //Cathode = General for xa=-(w+L)/4 and xb=(w+L)/4
      double f2_cath_lambda();
      double f2_cath_mu1();
      double f2_cath_mu2();
      double F2_cath_lambda(int n);
      double F2_cath_mu1(int n);
      double F2_cath_mu2(int n);
      double lambda_cath();
      double mu1_cath();
      double mu2_cath();
      
      //Top side
      double f2_top_mu2(double xa, double xb);
      double F2_top_mu2(int n, double xa, double xb);
      double mu2_top(double xa, double xb); 
      
      double t_f_top_mu2_0();
      double t_f_top_mu2_1();
      double t_f_top_mu2_2();
      double t_F_top_mu2_0(int n);
      double t_F_top_mu2_1(int n);
      double t_F_top_mu2_2(int n);
      
      double H_top_0();
      double H_top_1();
      double H_top_2();
      

    private:
      
      double d, p, w, L, g1, g2;
      int N;
      
      double _pi = 3.14159265359;      
      
      void flip_g1g2(); 
};

double GEM_Efficiencies_v3::Get_Coll_Driftions_g(double eta1, double eta2, double x1, double x2)
{
    double prefac = -1.0/(2.0*_pi*1.0);
    
    flip_g1g2();  
    double flux_anode = prefac*mu1_cath_g(-(w+L)/4.0,(w+L)/4.0)*eta2+prefac*lambda_cath_g(-(w+L)/4.0,(w+L)/4.0)*1.0+prefac*mu2_cath_g(-(w+L)/4.0,(w+L)/4.0)*eta1;
    double flux_D_tilde_ss = prefac*mu1_cath_g(0.0,x1)*eta2+prefac*lambda_cath_g(0.0,x1)*1.0+prefac*mu2_cath_g(0.0,x1)*eta1;
    double flux_A_ss = prefac*mu1_cath_g(x2,(w+L)/4.0)*eta2+prefac*lambda_cath_g(x2,(w+L)/4.0)*1.0+prefac*mu2_cath_g(x2,(w+L)/4.0)*eta1;
    flip_g1g2();
    
    double xb_bot = Get_xb_bot(eta1, eta2);
    double flux_A = prefac*(eta1-1.0)*mu2_top(-(w+L)/4.0, xb_bot)+eta2*(xb_bot+(w+L)/4.0);
    
    double coll_ion_drift = (0.5*flux_anode-flux_A-flux_D_tilde_ss)/(0.5*flux_anode-flux_A_ss-flux_D_tilde_ss);
   
    /*
    if ( coll_ion_drift <= 1.0 )
    {
	return coll_ion_drift;
    }else{
	return 1.0;
    }
    
    if ( coll_ion_drift >= 1.0 ){
	return 1.0;
    }else if ( coll_ion_drift >= 0.0 && coll_ion_drift < 1.0 ) {
	return coll_ion_drift;
    }else{
	return 0.0;
    }*/
    
    if ( coll_ion_drift < 0.0 )
    {
	return 0.0;
    }else{
	return coll_ion_drift;
    }
}


void GEM_Efficiencies_v3::flip_g1g2()
{
    double g1_old = g1;  
    g1 = g2;
    g2 = g1_old;
}

double GEM_Efficiencies_v3::Get_CollTop(double eta1, double eta2)
{
    return 2.0*_pi*(GetParameter_C7Bar(eta1,eta2)+GetParameter_C9Bar(eta1,eta2)*eta1+GetParameter_C8Bar(eta1,eta2)*eta2)/(GetParameter_C1()+GetParameter_C3()*eta1+GetParameter_C2()*eta2);
}

double GEM_Efficiencies_v3::Get_CollBot(double eta1, double eta2)
{
    return 2.0*_pi*(GetParameter_C7(eta1,eta2)+GetParameter_C8(eta1,eta2)*eta1+GetParameter_C9(eta1,eta2)*eta2)/(GetParameter_C1Bar()+GetParameter_C2Bar()*eta1+GetParameter_C3Bar()*eta2);
}

double GEM_Efficiencies_v3::Get_ExtrBot(double eta1, double eta2)
{
    return 2.0*_pi*(GetParameter_C7(eta1,eta2)+GetParameter_C8(eta1,eta2)*eta1+GetParameter_C9(eta1,eta2)*eta2)/(GetParameter_C4()+GetParameter_C5()*eta1+GetParameter_C6()*eta2);    
}

double GEM_Efficiencies_v3::Get_ExtrTop(double eta1, double eta2)
{
    return 2.0*_pi*(GetParameter_C7Bar(eta1,eta2)+GetParameter_C9Bar(eta1,eta2)*eta1+GetParameter_C8Bar(eta1,eta2)*eta2)/(GetParameter_C4()+GetParameter_C6()*eta1+GetParameter_C5()*eta2);
}

double GEM_Efficiencies_v3::Get_eta1k1(double eta2)
{
    return -1.0/(2.0*_pi)*H_top_0()*(1-eta2);
}

double GEM_Efficiencies_v3::Get_eta1k2(double eta2)
{
    return 1.0/(2.0*_pi)*((L-p)*(L-p)/4.0*H_top_2()+H_top_0())*(eta2-1.0);
}


double GEM_Efficiencies_v3::Get_eta2k1(double eta1)
{
    return -1.0/(2.0*_pi)*H_top_0()*(1-eta1);
}

double GEM_Efficiencies_v3::Get_eta2k2(double eta1)
{
    return 1.0/(2.0*_pi)*((L-p)*(L-p)/4.0*H_top_2()+H_top_0())*(eta1-1.0);
}

double GEM_Efficiencies_v3::Get_xb_top(double eta1, double eta2)
{    
    double result;
    
    if ( eta1 <= Get_eta1k1(eta2) )
    {
	result = -(L+w)/4.0;
    } 
    else if ( Get_eta1k1(eta2) < eta1 && eta1 < Get_eta1k2(eta2) )
    {
	result = -p/2.0 + sqrt(H_top_2()*(eta2-1.0)*(2.0*_pi*eta1+H_top_0()*(1.0-eta2)))/(H_top_2()*(eta2-1.0));
    }
    else
    {
	result = -L/2.0;
    }
    
    return result;
}

double GEM_Efficiencies_v3::Get_xb_bot(double eta1, double eta2)
{    
    double result;
    
    if ( eta2 <= Get_eta2k1(eta1) )
    {
	result = -(L+w)/4.0;
    } 
    else if ( Get_eta2k1(eta1) < eta2 && eta2 < Get_eta2k2(eta1) )
    {
	result = -p/2.0 + sqrt(H_top_2()*(eta1-1.0)*(2.0*_pi*eta2+H_top_0()*(1.0-eta1)))/(H_top_2()*(eta1-1.0));
    }
    else
    {
	result = -L/2.0;
    }
    
    return result;
}

double GEM_Efficiencies_v3::H_top_0()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++)
    {
	result += t_F_top_mu2_0(n);
    }
    
    result += t_f_top_mu2_0();
    
    return result;
}

double GEM_Efficiencies_v3::H_top_1()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += t_F_top_mu2_1(n);
    }
    
    result += t_f_top_mu2_1();
    
    return result;
}

double GEM_Efficiencies_v3::H_top_2()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += t_F_top_mu2_2(n);
    }
    
    result += t_f_top_mu2_2();
    
    return result;
}

//double GEM_Efficiencies_v3::t_F_top_mu2_0(int n, double d, double w, double L)
double GEM_Efficiencies_v3::t_F_top_mu2_0(int n)
{
    return atan(( ((n - 2) * L) +  ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1) / d / 0.2e1) + atan(( ((n - 2) * L) +  ((n - 1) * w) +  w 
    / 0.2e1 +  L / 0.2e1) / d / 0.2e1) - atan(( (n * L) +  ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1) / d / 0.2e1) - atan(( (n * L) +  
    ((n - 1) * w) +  w / 0.2e1 +  L / 0.2e1) / d / 0.2e1);
}

//double GEM_Efficiencies_v3::t_F_top_mu2_1(int n, double d, double w, double L)
double GEM_Efficiencies_v3::t_F_top_mu2_1(int n)
{
    return 0.1e1 / d / (pow( ((n - 2) * L) +  ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) - 0.1e1 / d 
    / (pow( ((n - 2) * L) +  ((n - 1) * w) +  w / 0.2e1 +  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) - 0.1e1 / d / (pow( 
    (n * L) +  ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) + 0.1e1 / d / (pow( (n * L) +  ((n - 1) * w) 
    +  w / 0.2e1 +  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1);
}

//double GEM_Efficiencies_v3::t_F_top_mu2_2(int n, double d, double w, double L)
double GEM_Efficiencies_v3::t_F_top_mu2_2(int n)
{
    return - (8 * (n - 2) * L + 8 * (n - 1) * w - 4 * w - 4 * L) / d / ( (4 *  pow( L,  2) *  pow( n,  2)) +  (8 * L *  pow( n,  2) * w) +  
    (4 *  pow( n,  2) *  pow( w,  2)) -  (20 *  pow( L,  2) * n) -  (32 * n * L * w) -  (12 * n *  pow( w,  2)) +  (25 *  pow( L,  2)) 
    +  (30 * L * w) + 0.16e2 * pow(d, 0.2e1) +  (9 *  pow( w,  2))) / (pow( ((n - 2) * L) +  ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1, 0.2e1) 
    * pow(d, -0.2e1) / 0.4e1 + 0.1e1) -  (8 * (n - 2) * L + 8 * (n - 1) * w + 4 * w + 4 * L) / d / ( (4 *  pow( L,  2) *  pow( n,  2)) +  
    (8 * L *  pow( n,  2) * w) +  (4 *  pow( n,  2) *  pow( w,  2)) -  (12 *  pow( L,  2) * n) -  (16 * n * L * w) -  (4 * n *  pow( w,  2)) 
    +  (9 *  pow( L,  2)) +  (6 * L * w) + 0.16e2 * pow(d, 0.2e1) +   pow( w,  2)) / (pow( ((n - 2) * L) +  ((n - 1) * w) +  w / 0.2e1 +  
    L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) +  (8 * n * L + 8 * (n - 1) * w - 4 * w - 4 * L) / d / ( (4 *  pow( L,  2) *  
    pow( n,  2)) +  (8 * L *  pow( n,  2) * w) +  (4 *  pow( n,  2) *  pow( w,  2)) -  (4 *  pow( L,  2) * n) -  (16 * n * L * w) -  
    (12 * n *  pow( w,  2)) +   pow( L,  2) +  (6 * L * w) + 0.16e2 * pow(d, 0.2e1) +  (9 *  pow( w,  2))) / (pow( (n * L) +  
    ((n - 1) * w) -  w / 0.2e1 -  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) +  (8 * n * L + 8 * (n - 1) * w + 4 * w + 4 * L) / d / 
    ( (4 *  pow( L,  2) *  pow( n,  2)) +  (8 * L *  pow( n,  2) * w) +  (4 *  pow( n,  2) *  pow( w,  2)) +  (4 *  pow( L,  2) * n) -  
    (4 * n *  pow( w,  2)) +   pow( L,  2) -  (2 * L * w) + 0.16e2 * pow(d, 0.2e1) +   pow( w,  2)) / (pow( (n * L) +  ((n - 1) * w) +  
    w / 0.2e1 +  L / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1);
}

//GEM_Efficiencies_v3::t_f_top_mu2_0(double d, double w, double L)
double GEM_Efficiencies_v3::t_f_top_mu2_0()
{
    return -atan((0.3e1 / 0.2e1 * L + w / 0.2e1) / d / 0.2e1) - atan((L / 0.2e1 - w / 0.2e1) / d / 0.2e1);
}

//GEM_Efficiencies_v3::t_f_top_mu2_0(double d, double w, double L)
double GEM_Efficiencies_v3::t_f_top_mu2_1()
{
    return 0.1e1 / d / (pow(0.3e1 / 0.2e1 * L + w / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) - 0.1e1 / d / 
    (pow(L / 0.2e1 - w / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1);
}
  
//GEM_Efficiencies_v3::t_f_top_mu2_0(double d, double w, double L)
double GEM_Efficiencies_v3::t_f_top_mu2_2()
{
    return (12 * L + 4 * w) / d / ( (9 *  pow( L,  2)) +  (6 * L * w) + 0.16e2 * pow(d, 0.2e1) +   pow( w,  2)) / (pow(0.3e1 / 0.2e1 
    *  L +  w / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1) +  (4 * L - 4 * w) / d / (  pow( L,  2) -  (2 * L * w) + 0.16e2 
    * pow(d, 0.2e1) +   pow( w,  2)) / (pow( L / 0.2e1 -  w / 0.2e1, 0.2e1) * pow(d, -0.2e1) / 0.4e1 + 0.1e1);
}

GEM_Efficiencies_v3::GEM_Efficiencies_v3(int _N)
{
    N = _N;
}

void GEM_Efficiencies_v3::SetGeometry(double _d, double _p, double _L, double _g1, double _g2)
{
    g1 = _g1;
    g2 = _g2;
    d = _d;
    p = _p;
    L = _L;
    w = 2*_p-_L;
}

void GEM_Efficiencies_v3::SetN(int _N)
{
    N = _N;
}

double GEM_Efficiencies_v3::GetParameter_C1()
{
    //return -lambda_cath();
    return -lambda_cath();
}

double GEM_Efficiencies_v3::GetParameter_C1Bar()
{
    flip_g1g2();
    return GetParameter_C1();
    flip_g1g2();
}

double GEM_Efficiencies_v3::GetParameter_C2()
{
    return -mu2_cath();
}

double GEM_Efficiencies_v3::GetParameter_C2Bar()
{
    flip_g1g2();
    return GetParameter_C2();
    flip_g1g2();
}

double GEM_Efficiencies_v3::GetParameter_C3()
{  
    return -mu1_cath();
}

double GEM_Efficiencies_v3::GetParameter_C3Bar()
{
    flip_g1g2();
    return GetParameter_C3();
    flip_g1g2();
}

double GEM_Efficiencies_v3::GetParameter_C4()
{
    return (mu2_top(-L/2.0, L/2.0) + _pi * L);
}

double GEM_Efficiencies_v3::GetParameter_C5()
{
    return -mu2_top(-L/2.0, L/2.0);
}

double GEM_Efficiencies_v3::GetParameter_C6()
{   
    return _pi*L;
}


double GEM_Efficiencies_v3::GetParameter_C7Bar_fromx(double xa, double xb)
{
    return (-1.0/(2.0*_pi))*(lambda_cath()+2.0*mu2_top(xa,xb));
    
}

double GEM_Efficiencies_v3::GetParameter_C8Bar_fromx(double xa, double xb)
{
    return (-1.0/(2.0*_pi))*(mu2_cath()-2.0*mu2_top(xa,xb));
}

double GEM_Efficiencies_v3::GetParameter_C9Bar_fromx(double xa, double xb)
{
    return (-1.0/(2.0*_pi))*(mu1_cath()+4*_pi*(xb-xa));
}

double GEM_Efficiencies_v3::GetParameter_C7Bar(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_top(eta1, eta2);
    
    return GetParameter_C7Bar_fromx(xa,xb);
}

double GEM_Efficiencies_v3::GetParameter_C7(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_bot(eta1, eta2);
    
    //flip g1<->g2 for calculation since c7bar(g1->g2) = c7
    flip_g1g2();
    double result = GetParameter_C7Bar_fromx(xa,xb);
    flip_g1g2();
    
    return result; 
}


double GEM_Efficiencies_v3::GetParameter_C8Bar(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_top(eta1, eta2);
    
    return GetParameter_C8Bar_fromx(xa,xb);
}

double GEM_Efficiencies_v3::GetParameter_C8(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_bot(eta1, eta2);
    
    flip_g1g2();
    double result = GetParameter_C8Bar_fromx(xa,xb);
    flip_g1g2();
    
    return result;
}

double GEM_Efficiencies_v3::GetParameter_C9Bar(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_top(eta1, eta2);
        
    return GetParameter_C9Bar_fromx(xa,xb);
}

double GEM_Efficiencies_v3::GetParameter_C9(double eta1, double eta2)
{
    double xa = -(w+L)/4.0;
    double xb = Get_xb_bot(eta1, eta2);
    
    flip_g1g2();
    double result = GetParameter_C9Bar_fromx(xa,xb);
    flip_g1g2();
    
    return result;
}

double GEM_Efficiencies_v3::mu1_cath()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += F2_cath_mu1(n);
    }
    
    result += f2_cath_mu1();
    
    return result;
}

double GEM_Efficiencies_v3::mu2_cath()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += F2_cath_mu2(n);
    }
    
    result += f2_cath_mu2();
    
    return result;
}

double GEM_Efficiencies_v3::lambda_cath()
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += F2_cath_lambda(n);
    }
    
    result += f2_cath_lambda();
    
    return result;
}

double GEM_Efficiencies_v3::mu2_top(double xa, double xb)
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += F2_top_mu2(n, xa, xb);
    }
    
    result += f2_top_mu2(xa, xb);
    
    return result;
}

//double GEM_Efficiencies_v3::f2_top_mu2(double d, double L, double xa, double xb)
double GEM_Efficiencies_v3::f2_top_mu2(double xa, double xb)
{
    return atan( ((L - 2 * xa) / d) / 0.2e1) *  xa - atan( ((L - 2 * xa) / d) / 0.2e1) *  L / 0.2e1 +  d * log( ( pow( L,  2) - 4 * L * xa + 
    4 *  pow( d,  2) + 4 * xa * xa)) / 0.2e1 + atan( ((L + 2 * xa) / d) / 0.2e1) *  xa + atan( ((L + 2 * xa) / d) / 0.2e1) *  L / 0.2e1 -  d * log( ( pow( L,  2) + 
    4 * L * xa + 4 *  pow( d,  2) + 4 * xa * xa)) / 0.2e1 - atan( ((-2 * xb + L) / d) / 0.2e1) *  xb + atan( ((-2 * xb + L) / d) / 0.2e1) *  L / 0.2e1 -  d * 
    log( ( pow( L,  2) - 4 * L * xb + 4 *  pow( d,  2) + 4 * xb * xb)) / 0.2e1 - atan( ((2 * xb + L) / d) / 0.2e1) *  xb - atan( ((2 * xb + L) / d) / 0.2e1) *  
    L / 0.2e1 +  d * log( ( pow( L,  2) + 4 * L * xb + 4 *  pow( d,  2) + 4 * xb * xb)) / 0.2e1;
}

//double GEM_Efficiencies_v3::F2_top_mu2(int n, double d, double w, double L, double xa, double xb)
double GEM_Efficiencies_v3::F2_top_mu2(int n, double xa, double xb)
{
    return atan( ((n * L + n * w - w + 2 * xa) / d) / 0.2e1) *  xa - atan( ((n * L + n * w - w + 2 * xa) / d) / 0.2e1) *  w / 0.2e1 + 
    atan( ((n * L + n * w - w - 2 * xa) / d) / 0.2e1) *  xa + atan( ((n * L + n * w - w - 2 * xa) / d) / 0.2e1) *  w / 0.2e1 + atan( ((n * L + n * w - 2 * L - 
    w + 2 * xb) / d) / 0.2e1) *  xb - atan( ((n * L + n * w - 2 * L - w + 2 * xb) / d) / 0.2e1) *  L - atan( ((n * L + n * w - 2 * L - w + 2 * xb) / d) / 0.2e1) *  
    w / 0.2e1 + atan( ((n * L + n * w - 2 * L - w - 2 * xb) / d) / 0.2e1) *  xb + atan( ((n * L + n * w - 2 * L - w - 2 * xb) / d) / 0.2e1) *  L + atan( ((n * L + n * 
    w - 2 * L - w - 2 * xb) / d) / 0.2e1) *  w / 0.2e1 - atan( ((n * L + n * w - w + 2 * xb) / d) / 0.2e1) *  xb + atan( ((n * L + n * w - w + 2 * xb) / d) / 0.2e1) *  
    w / 0.2e1 - atan( ((n * L + n * w - w - 2 * xb) / d) / 0.2e1) *  xb - atan( ((n * L + n * w - w - 2 * xb) / d) / 0.2e1) *  w / 0.2e1 - atan( ((n * L + n * w - 2 * 
    L - w + 2 * xa) / d) / 0.2e1) *  xa + atan( ((n * L + n * w - 2 * L - w + 2 * xa) / d) / 0.2e1) *  L + atan( ((n * L + n * w - 2 * L - w + 2 * xa) / d) / 0.2e1) *  
    w / 0.2e1 - atan( ((n * L + n * w - 2 * L - w - 2 * xa) / d) / 0.2e1) *  xa - atan( ((n * L + n * w - 2 * L - w - 2 * xa) / d) / 0.2e1) *  L - atan( ((n * L + n * w 
    - 2 * L - w - 2 * xa) / d) / 0.2e1) *  w / 0.2e1 + atan( ((n * L + n * w - 2 * L - w - 2 * xa) / d) / 0.2e1) *  n *  L / 0.2e1 + atan( ((n * L + n * w - 2 * L - w - 
    2 * xa) / d) / 0.2e1) *  n *  w / 0.2e1 + atan( ((n * L + n * w - w + 2 * xa) / d) / 0.2e1) *  n *  w / 0.2e1 + atan( ((n * L + n * w - w + 2 * xa) / d) / 0.2e1) *  
    n *  L / 0.2e1 - atan( ((n * L + n * w - w - 2 * xa) / d) / 0.2e1) *  n *  L / 0.2e1 - atan( ((n * L + n * w - w - 2 * xa) / d) / 0.2e1) *  n *  w / 0.2e1 + atan( 
    ((n * L + n * w - 2 * L - w + 2 * xb) / d) / 0.2e1) *  n *  L / 0.2e1 + atan( ((n * L + n * w - 2 * L - w + 2 * xb) / d) / 0.2e1) *  n *  w / 0.2e1 - atan( ((n * L + 
    n * w - 2 * L - w - 2 * xb) / d) / 0.2e1) *  n *  L / 0.2e1 - atan( ((n * L + n * w - 2 * L - w - 2 * xb) / d) / 0.2e1) *  n *  w / 0.2e1 - atan( ((n * L + n * w - w + 
    2 * xb) / d) / 0.2e1) *  n *  L / 0.2e1 - atan( ((n * L + n * w - w + 2 * xb) / d) / 0.2e1) *  n *  w / 0.2e1 + atan( ((n * L + n * w - w - 2 * xb) / d) / 0.2e1) *  n 
    *  L / 0.2e1 + atan( ((n * L + n * w - w - 2 * xb) / d) / 0.2e1) *  n *  w / 0.2e1 - atan( ((n * L + n * w - 2 * L - w + 2 * xa) / d) / 0.2e1) *  n *  L / 0.2e1 - 
    atan( ((n * L + n * w - 2 * L - w + 2 * xa) / d) / 0.2e1) *  n *  w / 0.2e1 +  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  
    pow( w,  2) - 4 *  pow( L,  2) * n - 6 * n * L * w - 4 * n * L * xb - 2 * n *  pow( w,  2) - 4 * xb * w * n + 4 *  pow( L,  2) + 4 * L * w + 8 * L * xb + 4 *  
    pow( d,  2) +  pow( w,  2) + 4 * xb * w + 4 * xb * xb)) / 0.2e1 -  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 
    2 * n * L * w + 4 * n * L * xa - 2 * n *  pow( w,  2) + 4 * w * xa * n + 4 *  pow( d,  2) +  pow( w,  2) - 4 * w * xa + 4 * xa * xa)) / 0.2e1 -  d * log( ( pow( L,  2) 
    *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 2 * n * L * w - 4 * n * L * xb - 2 * n *  pow( w,  2) - 4 * xb * w * n + 4 *  pow( d,  2) +  
    pow( w,  2) + 4 * xb * w + 4 * xb * xb)) / 0.2e1 -  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) 
    * n - 6 * n * L * w - 4 * n * L * xa - 2 * n *  pow( w,  2) - 4 * w * xa * n + 4 *  pow( L,  2) + 4 * L * w + 8 * L * xa + 4 *  pow( d,  2) +  pow( w,  2) + 4 * w * 
    xa + 4 * xa * xa)) / 0.2e1 +  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) * n - 6 * n * L * w + 
    4 * n * L * xa - 2 * n *  pow( w,  2) + 4 * w * xa * n + 4 *  pow( L,  2) + 4 * L * w - 8 * L * xa + 4 *  pow( d,  2) +  pow( w,  2) - 4 * w * xa + 4 * xa * xa)) / 0.2e1 
    -  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) * n - 6 * n * L * w + 4 * n * L * xb - 2 * n *  
    pow( w,  2) + 4 * xb * w * n + 4 *  pow( L,  2) + 4 * L * w - 8 * L * xb + 4 *  pow( d,  2) +  pow( w,  2) - 4 * xb * w + 4 * xb * xb)) / 0.2e1 +  d * log( 
    ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) *  pow( w,  2) - 2 * n * L * w - 4 * n * L * xa - 2 * n *  pow( w,  2) - 4 * w * xa * 
    n + 4 *  pow( d,  2) +  pow( w,  2) + 4 * w * xa + 4 * xa * xa)) / 0.2e1 +  d * log( ( pow( L,  2) *  pow( n,  2) + 2 * L *  pow( n,  2) * w +  pow( n,  2) 
    *  pow( w,  2) - 2 * n * L * w + 4 * n * L * xb - 2 * n *  pow( w,  2) + 4 * xb * w * n + 4 *  pow( d,  2) +  pow( w,  2) - 4 * xb * w + 4 * xb * xb)) / 0.2e1;
    
}

//double GEM_Efficiencies_v3::f2_cath_mu2(double d, double w, double L, double g1)
double GEM_Efficiencies_v3::f2_cath_mu2()
{
    return -0.3e1 / 0.2e1 *  L * atan( ((3 * L + w) / (d + g1)) / 0.2e1) +  L * atan( ((-w + L) / (d + g1)) / 0.2e1) / 0.2e1 - 
    atan( ((3 * L + w) / (d + g1)) / 0.2e1) *  w / 0.2e1 +  d * log( (9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) 
    +  pow( w,  2))) / 0.2e1 +  g1 * log( (9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 - 
    atan( ((-w + L) / (d + g1)) / 0.2e1) *  w / 0.2e1 -  d * log( ( pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  
    pow( w,  2))) / 0.2e1 -  g1 * log( ( pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1;
}

//double GEM_Efficiencies_v3::F2_cath_mu2(int n, double d, double w, double L, double g1)
double GEM_Efficiencies_v3::F2_cath_mu2(int n)
{
    return - L * atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (d + g1)) / 0.2e1) *  n +  L * atan( ((2 * n * L + 2 * n * w - L - 3 * w) 
    / (d + g1)) / 0.2e1) *  n - atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (d + g1)) / 0.2e1) *  n *  w + atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (d + g1)) 
    / 0.2e1) *  n *  w - atan( ((2 * n * L + 2 * n * w + L - w) / (d + g1)) / 0.2e1) *  n *  w + atan( ((2 * n * L + 2 * n * w - L - 3 * w) / (d + g1)) / 0.2e1) *  n 
    *  w +  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) + 4 *  pow( L,  2) * n - 4 * n *  pow( w,  2) +  
    pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 +  g1 * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  
    pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) + 4 *  pow( L,  2) * n - 4 * n *  pow( w,  2) +  pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 
    *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 -  L * atan( ((2 * n * L + 2 * n * w + L - w) / (d + g1)) / 0.2e1) *  n -  L * atan( ((2 * n * L + 2 * n * w + L - w) / 
    (d + g1)) / 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (d + g1)) / 0.2e1) *  w - atan( ((2 * n * L + 2 * n * w - 3 * L - w) / 
    (d + g1)) / 0.2e1) *  w / 0.2e1 + atan( ((2 * n * L + 2 * n * w + L - w) / (d + g1)) / 0.2e1) *  w / 0.2e1 - 0.3e1 / 0.2e1 * atan( ((2 * n * L + 2 * n * w - L - 3 * w) / 
    (d + g1)) / 0.2e1) *  w +  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 20 *  pow( L,  2) * n - 32 * n * 
    L * w - 12 * n *  pow( w,  2) + 25 *  pow( L,  2) + 30 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 +  g1 * log( (4 *  
    pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 20 *  pow( L,  2) * n - 32 * n * L * w - 12 * n *  pow( w,  2) + 25 *  
    pow( L,  2) + 30 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 -  L * atan( ((2 * n * L + 2 * n * w - L - 3 * w) / 
    (d + g1)) / 0.2e1) / 0.2e1 + 0.5e1 / 0.2e1 *  L * atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (d + g1)) / 0.2e1) -  g1 * log( (4 *  pow( L,  2) *  pow( n,  2) + 
    8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) * n - 16 * n * L * w - 12 * n *  pow( w,  2) +  pow( L,  2) + 6 * L * w + 4 *  
    pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 +  L * atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (d + g1)) / 0.2e1) *  n - 0.3e1 / 0.2e1 
    *  L * atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (d + g1)) / 0.2e1) -  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  
    pow( w,  2) - 4 *  pow( L,  2) * n - 16 * n * L * w - 12 * n *  pow( w,  2) +  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  
    pow( w,  2))) / 0.2e1 -  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 12 *  pow( L,  2) * n - 16 * n * 
    L * w - 4 * n *  pow( w,  2) + 9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 -  g1 * log( (4 *  
    pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 12 *  pow( L,  2) * n - 16 * n * L * w - 4 * n *  pow( w,  2) + 9 
    *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) + 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1;
}

//double GEM_Efficiencies_v3::f2_cath_mu1(double d, double w, double L, double g1)
double GEM_Efficiencies_v3::f2_cath_mu1()
{
    return -0.3e1 / 0.2e1 *  L * atan( ((3 * L + w) / (-g1 + d)) / 0.2e1) +  L * atan( ((-w + L) / (-g1 + d)) / 0.2e1) / 0.2e1 
    -  L * 0.3141592654e1 - atan( ((3 * L + w) / (-g1 + d)) / 0.2e1) *  w / 0.2e1 +  d * log( (9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) - 
    8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 -  g1 * log( (9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  
    pow( g1,  2) +  pow( w,  2))) / 0.2e1 - atan( ((-w + L) / (-g1 + d)) / 0.2e1) *  w / 0.2e1 -  d * log( ( pow( L,  2) - 2 * L * w + 4 *  
    pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 +  g1 * log( ( pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) - 8 * g1 * 
    d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 - 0.3141592654e1 *  w;
}

//double GEM_Efficiencies_v3::F2_cath_mu1(int n, double d, double w, double L, double g1)
double GEM_Efficiencies_v3::F2_cath_mu1(int n)
{
    return -0.3e1 / 0.2e1 *  L * atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (-g1 + d)) / 0.2e1) +  L * 
    atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (-g1 + d)) / 0.2e1) *  n +  L * atan( ((2 * n * L + 2 * n * w - L - 3 * w) / (-g1 + d)) / 0.2e1) *  
    n - atan( ((2 * n * L + 2 * n * w + L - w) / (-g1 + d)) / 0.2e1) *  n *  w - atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (-g1 + d)) / 0.2e1) 
    *  n *  w + atan( ((2 * n * L + 2 * n * w - 3 * L - w) / (-g1 + d)) / 0.2e1) *  n *  w + atan( ((2 * n * L + 2 * n * w - L - 3 * w) / (-g1 + d)) / 0.2e1) 
    *  n *  w -  L * atan( ((2 * n * L + 2 * n * w + L - w) / (-g1 + d)) / 0.2e1) *  n + 0.5e1 / 0.2e1 *  L * atan( ((2 * n * L + 2 * n * w - 5 * L 
    - 3 * w) / (-g1 + d)) / 0.2e1) -  g1 * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 20 *  
    pow( L,  2) * n - 32 * n * L * w - 12 * n *  pow( w,  2) + 25 *  pow( L,  2) + 30 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) + 9 
    *  pow( w,  2))) / 0.2e1 +  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 20 *  pow( L,  2) 
    * n - 32 * n * L * w - 12 * n *  pow( w,  2) + 25 *  pow( L,  2) + 30 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 
    -  L * atan( ((2 * n * L + 2 * n * w - L - 3 * w) / (-g1 + d)) / 0.2e1) / 0.2e1 + atan( ((2 * n * L + 2 * n * w + L - w) / (-g1 + d)) / 0.2e1) *  w / 
    0.2e1 + 0.3e1 / 0.2e1 * atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (-g1 + d)) / 0.2e1) *  w - atan( ((2 * n * L + 2 * n * w - 3 * L - w) / 
    (-g1 + d)) / 0.2e1) *  w / 0.2e1 - 0.3e1 / 0.2e1 * atan( ((2 * n * L + 2 * n * w - L - 3 * w) / (-g1 + d)) / 0.2e1) *  w +  g1 * log( (4 *  pow( L,  2) 
    *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) * n - 16 * n * L * w - 12 * n *  pow( w,  2) +  
    pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 -  d * log( (4 *  pow( L,  2) *  pow( n,  2) 
    + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 4 *  pow( L,  2) * n - 16 * n * L * w - 12 * n *  pow( w,  2) +  pow( L,  2) + 6 * L * 
    w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) + 9 *  pow( w,  2))) / 0.2e1 -  L * atan( ((2 * n * L + 2 * n * w - 5 * L - 3 * w) / (-g1 + d)) / 0.2e1) 
    *  n -  L * atan( ((2 * n * L + 2 * n * w + L - w) / (-g1 + d)) / 0.2e1) / 0.2e1 -  g1 * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 
    *  pow( n,  2) *  pow( w,  2) + 4 *  pow( L,  2) * n - 4 * n *  pow( w,  2) +  pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) +  
    pow( w,  2))) / 0.2e1 +  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) + 4 *  pow( L,  2) * n - 4 * 
    n *  pow( w,  2) +  pow( L,  2) - 2 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 +  g1 * log( (4 *  pow( L,  2) *  
    pow( n,  2) + 8 * L *  pow( n,  2) * w + 4 *  pow( n,  2) *  pow( w,  2) - 12 *  pow( L,  2) * n - 16 * n * L * w - 4 * n *  pow( w,  2) + 9 *  pow( L,  2) 
    + 6 * L * w + 4 *  pow( d,  2) - 8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1 -  d * log( (4 *  pow( L,  2) *  pow( n,  2) + 8 * L *  pow( n,  2) 
    * w + 4 *  pow( n,  2) *  pow( w,  2) - 12 *  pow( L,  2) * n - 16 * n * L * w - 4 * n *  pow( w,  2) + 9 *  pow( L,  2) + 6 * L * w + 4 *  pow( d,  2) - 
    8 * g1 * d + 4 *  pow( g1,  2) +  pow( w,  2))) / 0.2e1;
}

//double GEM_Efficiencies_v3::f2_cath_lambda(double d, double w, double L, double g1)
double GEM_Efficiencies_v3::f2_cath_lambda()
{
    return -d * log(0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 
    pow(w, 0.2e1)) / 0.2e1 + g1 * log(0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 
    0.2e1 + atan((-w + L) / (-g1 + d) / 0.2e1) * w / 0.2e1 + 0.3e1 / 0.2e1 * L * atan((0.3e1 * L + w) / (d + g1) / 0.2e1) + 0.3e1 / 0.2e1 * L * atan((0.3e1 * L + w) / 
    (-g1 + d) / 0.2e1) - d * log(0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 
    0.2e1 - g1 * log(0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - L * atan((-w + L) / 
    (d + g1) / 0.2e1) / 0.2e1 + atan((0.3e1 * L + w) / (d + g1) / 0.2e1) * w / 0.2e1 - L * atan((-w + L) / (-g1 + d) / 0.2e1) / 0.2e1 + d * log(pow(L, 0.2e1) - 
    0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - g1 * log(pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) 
    - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 + atan((0.3e1 * L + w) / (-g1 + d) / 0.2e1) * w / 0.2e1 + atan((-w + L) / (d + g1) / 0.2e1) * w / 0.2e1 
    + d * log(pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 + g1 * log(pow(L, 0.2e1) - 0.2e1 
    * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1;
}

//double GEM_Efficiencies_v3::F2_cath_lambda(int n, double d, double w, double L, double g1)
double GEM_Efficiencies_v3::F2_cath_lambda(int n)
{
    return -d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) + 0.4e1 * 
    pow(L, 0.2e1) * n - 0.4e1 * n * pow(w, 0.2e1) + pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - 
    g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) + 0.4e1 * pow(L, 0.2e1) * n - 0.4e1 * n * pow(w, 0.2e1) 
    + pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - 0.3e1 / 0.2e1 * atan((0.2e1 * n * L + 0.2e1 * 
    n * w - 0.5e1 * L - 0.3e1 * w) / (d + g1) / 0.2e1) * w + atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / (d + g1) / 0.2e1) * w / 0.2e1 - 
    atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (d + g1) / 0.2e1) * w / 0.2e1 + 0.3e1 / 0.2e1 * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / (d + g1) / 0.2e1) * w - 
    d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.20e2 * pow(L, 0.2e1) * n - 0.32e2 * n * L * w - 
    0.12e2 * n * pow(w, 0.2e1) + 0.25e2 * pow(L, 0.2e1) + 0.30e2 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 - 
    g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.20e2 * pow(L, 0.2e1) * n - 0.32e2 * n * L * w - 
    0.12e2 * n * pow(w, 0.2e1) + 0.25e2 * pow(L, 0.2e1) + 0.30e2 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 + 
    g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 
    0.12e2 * n * pow(w, 0.2e1) + pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 + d * 
    log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 
    0.12e2 * n * pow(w, 0.2e1) + pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 + d * 
    log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.12e2 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w -
    0.4e1 * n * pow(w, 0.2e1) + 0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 + g1 * 
    log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.12e2 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 
    0.4e1 * n * pow(w, 0.2e1) + 0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) + 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 + L * 
    atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (-g1 + d) / 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / (-g1 + d) / 0.2e1) 
    - 0.5e1 / 0.2e1 * L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / (-g1 + d) / 0.2e1) + g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * 
    L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.20e2 * pow(L, 0.2e1) * n - 0.32e2 * n * L * w - 0.12e2 * n * pow(w, 0.2e1) + 0.25e2 * pow(L, 0.2e1) + 
    0.30e2 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 - d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 
    0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.20e2 * pow(L, 0.2e1) * n - 0.32e2 * n * L * w - 0.12e2 * n * pow(w, 0.2e1) + 0.25e2 * 
    pow(L, 0.2e1) + 0.30e2 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 - atan((0.2e1 * n * L + 0.2e1 * 
    n * w + L - w) / (-g1 + d) / 0.2e1) * w / 0.2e1 - 0.3e1 / 0.2e1 * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / (-g1 + d) / 0.2e1) * 
    w + atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / (-g1 + d) / 0.2e1) * w / 0.2e1 + 0.3e1 / 0.2e1 * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / 
    (-g1 + d) / 0.2e1) * w - g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 0.12e2 * n * pow(w, 0.2e1) + pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * 
    pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 + d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 0.12e2 * n * pow(w, 0.2e1) + pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * 
    g1 * d + 0.4e1 * pow(g1, 0.2e1) + 0.9e1 * pow(w, 0.2e1)) / 0.2e1 + g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * 
    pow(n, 0.2e1) * pow(w, 0.2e1) + 0.4e1 * pow(L, 0.2e1) * n - 0.4e1 * n * pow(w, 0.2e1) + pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * 
    d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * 
    pow(w, 0.2e1) + 0.4e1 * pow(L, 0.2e1) * n - 0.4e1 * n * pow(w, 0.2e1) + pow(L, 0.2e1) - 0.2e1 * L * w + 0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * 
    pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - g1 * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * w + 0.4e1 * pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.12e2 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 0.4e1 * n * pow(w, 0.2e1) + 0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 0.4e1 * 
    pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 + d * log(0.4e1 * pow(L, 0.2e1) * pow(n, 0.2e1) + 0.8e1 * L * pow(n, 0.2e1) * 
    w + 0.4e1 * pow(n, 0.2e1) * pow(w, 0.2e1) - 0.12e2 * pow(L, 0.2e1) * n - 0.16e2 * n * L * w - 0.4e1 * n * pow(w, 0.2e1) + 0.9e1 * pow(L, 0.2e1) + 0.6e1 * L * w + 
    0.4e1 * pow(d, 0.2e1) - 0.8e1 * g1 * d + 0.4e1 * pow(g1, 0.2e1) + pow(w, 0.2e1)) / 0.2e1 - 0.5e1 / 0.2e1 * L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / 
    (d + g1) / 0.2e1) + L * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / (d + g1) / 0.2e1) / 0.2e1 + L * atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (d + g1) 
    / 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / (d + g1) / 0.2e1) + L * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) 
    / (-g1 + d) / 0.2e1) / 0.2e1 + atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / (-g1 + d) / 0.2e1) * n * w - atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * 
    L - w) / (d + g1) / 0.2e1) * n * w + atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (-g1 + d) / 0.2e1) * n * w - L * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / 
    (-g1 + d) / 0.2e1) * n - L * atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / (d + g1) / 0.2e1) * n + L * atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / 
    (-g1 + d) / 0.2e1) * n + L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / (-g1 + d) / 0.2e1) * n - atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / 
    (-g1 + d) / 0.2e1) * n * w + atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / (d + g1) / 0.2e1) * n * w - atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / 
    (-g1 + d) / 0.2e1) * n * w - L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / (-g1 + d) / 0.2e1) * n - atan((0.2e1 * n * L + 0.2e1 * n * w - L - 0.3e1 * w) / 
    (d + g1) / 0.2e1) * n * w + L * atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (d + g1) / 0.2e1) * n + L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.5e1 * L - 0.3e1 * w) / 
    (d + g1) / 0.2e1) * n + atan((0.2e1 * n * L + 0.2e1 * n * w + L - w) / (d + g1) / 0.2e1) * n * w - L * atan((0.2e1 * n * L + 0.2e1 * n * w - 0.3e1 * L - w) / 
    (d + g1) / 0.2e1) * n;
}

double GEM_Efficiencies_v3::fg_cath_lambda(double xa, double xb)
{
    return L * atan( ((-2 * xa + L) / (d + g1))) / 0.2e1 -  L * atan( ((2 * xa + L) / (d + g1))) / 0.2e1 -  L * atan( ((-2 * xb + L) / (d + g1))) / 0.2e1 +  
    L * atan( ((2 * xb + L) / (d + g1))) / 0.2e1 - atan( ((-2 * xa + L) / (d + g1))) *  xa - atan( ((2 * xa + L) / (d + g1))) *  xa + atan( ((-2 * xb + L) / 
    (d + g1))) *  xb + atan( ((2 * xb + L) / (d + g1))) *  xb -  d * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) 
    / 0.4e1 -  g1 * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 +  d * log( ( pow( L,  2) + 4 * L * xa +  
    pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 +  g1 * log( ( pow( L,  2) + 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  
    pow( g1,  2) + 4 * xa * xa)) / 0.4e1 +  d * log( ( pow( L,  2) - 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 +  
    g1 * log( ( pow( L,  2) - 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 -  d * log( ( pow( L,  2) + 4 * L * xb +  
    pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 -  g1 * log( ( pow( L,  2) + 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 
    4 * xb * xb)) / 0.4e1 +  d * log( ( pow( L,  2) - 4 * L * xb +  pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 +  g1 * 
    log( ( pow( L,  2) + 4 * L * xb +  pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 -  d * log( ( pow( L,  2) + 4 * L * xb +  
    pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 -  d * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) - 2 * d * g1 +  
    pow( g1,  2) + 4 * xa * xa)) / 0.4e1 -  g1 * log( ( pow( L,  2) + 4 * L * xa +  pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 +  
    d * log( ( pow( L,  2) + 4 * L * xa +  pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 -  g1 * log( ( pow( L,  2) - 4 * L * xb +  
    pow( d,  2) - 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 +  g1 * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) - 2 * d * g1 +  
    pow( g1,  2) + 4 * xa * xa)) / 0.4e1 - atan( ((-2 * xa + L) / (-g1 + d))) *  xa - atan( ((2 * xa + L) / (-g1 + d))) *  xa + atan( ((-2 * xb + L) / 
    (-g1 + d))) *  xb + atan( ((2 * xb + L) / (-g1 + d))) *  xb +  L * atan( ((-2 * xa + L) / (-g1 + d))) / 0.2e1 -  L * atan( ((2 * xa + L) / 
    (-g1 + d))) / 0.2e1 -  L * atan( ((-2 * xb + L) / (-g1 + d))) / 0.2e1 +  L * atan( ((2 * xb + L) / (-g1 + d))) / 0.2e1;
}

double GEM_Efficiencies_v3::fg_cath_mu1(double xa, double xb)
{
    return 0.2e1 * 0.3141592654e1 * xa - 0.2e1 * 0.3141592654e1 * xb - d * log(pow(L, 0.2e1) - 0.4e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) + 0.4e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) + 0.4e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) - 0.4e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) + 0.4e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) + 0.4e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) - 0.4e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) - 0.4e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 0.4e1 * xa * xa) / 0.4e1 + atan((-0.2e1 * xa + L) / (-g1 + d)) * xa + atan((0.2e1 * xa + L) / (-g1 + d)) * xa - atan((-0.2e1 * xb + L) / (-g1 + d)) * xb - atan((0.2e1 * xb + L) / (-g1 + d)) * xb - L * atan((-0.2e1 * xa + L) / (-g1 + d)) / 0.2e1 + L * atan((0.2e1 * xa + L) / (-g1 + d)) / 0.2e1 + L * atan((-0.2e1 * xb + L) / (-g1 + d)) / 0.2e1 - L * atan((0.2e1 * xb + L) / (-g1 + d)) / 0.2e1;
}

double GEM_Efficiencies_v3::fg_cath_mu2(double xa, double xb)
{
    return - L * atan( ((-2 * xa + L) / (d + g1))) / 0.2e1 +  L * atan( ((2 * xa + L) / (d + g1))) / 0.2e1 +  L * atan( ((-2 * xb + L) / (d + g1))) / 0.2e1 -  L * atan( ((2 * xb + L) / (d + g1))) / 0.2e1 + atan( ((-2 * xa + L) / (d + g1))) *  xa +  d * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 +  g1 * log( ( pow( L,  2) - 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 + atan( ((2 * xa + L) / (d + g1))) *  xa -  d * log( ( pow( L,  2) + 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 -  g1 * log( ( pow( L,  2) + 4 * L * xa +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xa * xa)) / 0.4e1 - atan( ((-2 * xb + L) / (d + g1))) *  xb -  d * log( ( pow( L,  2) - 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 -  g1 * log( ( pow( L,  2) - 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 - atan( ((2 * xb + L) / (d + g1))) *  xb +  d * log( ( pow( L,  2) + 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1 +  g1 * log( ( pow( L,  2) + 4 * L * xb +  pow( d,  2) + 2 * d * g1 +  pow( g1,  2) + 4 * xb * xb)) / 0.4e1;
}

double GEM_Efficiencies_v3::Fg_cath_lambda(int n, double xa, double xb)
{
    return -d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * 
    L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 
    0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * 
    w + 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 
    0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 0.4e1 * 
    pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 
    0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * 
    L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) - 
    0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * 
    pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 
    pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 
    0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * 
    xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 
    0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 
    0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 
    0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 0.4e1 * 
    pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 
    0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * 
    xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * 
    xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 
    0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * 
    w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 
    0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 
    0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + 
    pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * 
    w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) - 
    0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / 
    (-g1 + d)) * n / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (-g1 + d)) - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) + L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) + 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * xb + 
    atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * xb - atan((n * L + w * n - w + 
    0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * xb - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (-g1 + d)) * xa - atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * 
    w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * xa + atan((n * L + w * n - w + 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 - atan((n * L + 
    w * n - w + 0.2e1 * xa) / (-g1 + d)) * xa - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * 
    L - w - 0.2e1 * xb) / (-g1 + d)) * xb + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - 
    w - 0.2e1 * xa) / (-g1 + d)) * xa + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 
    0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 
    0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + 
    pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * 
    w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) + 
    0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * 
    pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 
    pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * 
    xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + 
    pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 
    0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * 
    xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 
    0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * 
    xa + 0.4e1 * xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * 
    n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * 
    w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * 
    L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 
    0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * 
    L * w - 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * 
    log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 
    0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * 
    d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * 
    w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * 
    xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * 
    xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * 
    n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + 
    pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 
    0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * 
    w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * 
    pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 
    0.4e1 - L * atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * n / 0.2e1 + L * atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * n / 
    0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 + L * atan((n * L + w * n - w - 0.2e1 * xa) / 
    (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - w + 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w - 
    0.2e1 * xb) / (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 - 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xb) / 
    (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w + 
    0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1 - atan((n * L + w * n - w + 
    0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 - atan((n * L + 
    w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) * n / 0.2e1 - 
    L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) * n / 0.2e1 + L * atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * n / 
    0.2e1 - L * atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * n / 0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * 
    n / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * n / 0.2e1 - L * atan((n * L + w * n - w - 0.2e1 * xb) / 
    (d + g1)) * n / 0.2e1 + L * atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * n / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / 
    (d + g1)) * n * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) * n * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xa) / 
    (d + g1)) * n * w / 0.2e1 - atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / 
    (d + g1)) * n * w / 0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / 
    (d + g1)) + L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (d + g1)) - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / 
    (d + g1)) * xb + atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / 
    (d + g1)) * xb - atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xa) / 
    (d + g1)) * xa + atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - w + 0.2e1 * xa) / 
    (d + g1)) * xa + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - 
    w - 0.2e1 * xa) / (d + g1)) * xa - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * 
    L - w + 0.2e1 * xa) / (d + g1)) * xa + atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xb) / 
    (d + g1)) * xb - atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * xb + 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * n * w / 
    0.2e1 + atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * n * w / 0.2e1;
}

double GEM_Efficiencies_v3::Fg_cath_mu1(int n, double xa, double xb)
{
    return d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 
    0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 
    0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 
    0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 
    0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * 
    L * w + 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + d * 
    log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * 
    w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + 
    pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * 
    n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * 
    d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * 
    pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 
    pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * 
    pow(w, 0.2e1) - 0.4e1 * n * w * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 
    0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * 
    L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * 
    xb + 0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.2e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * 
    w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * 
    w * xa + 0.4e1 * xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 0.4e1 * pow(L, 0.2e1) + 
    0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + 
    g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * 
    xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * 
    xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.2e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * 
    n * w * xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 
    0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 
    0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * 
    xb + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) - 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * n / 
    0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (-g1 + d)) + L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / 
    (-g1 + d)) - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / 
    (-g1 + d)) * xb - atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * xb + 
    atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 - atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * xb + 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (-g1 + d)) * xa + atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xa) / 
    (-g1 + d)) * xa - atan((n * L + w * n - w + 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xa) / 
    (-g1 + d)) * xa + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * 
    L - w - 0.2e1 * xb) / (-g1 + d)) * xb - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * w / 0.2e1 - 
    atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * xa + L * atan((n * L + w * n - w - 0.2e1 * xb) / 
    (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - 
    0.2e1 * L - w + 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 + 
    L * atan((n * L + w * n - w + 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / 
    (-g1 + d)) * n / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * n / 0.2e1 + atan((n * L + w * n - 
    0.2e1 * L - w + 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 - 
    atan((n * L + w * n - w + 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (-g1 + d)) * 
    n * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xa) / 
    (-g1 + d)) * n * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (-g1 + d)) * n * w / 0.2e1 + atan((n * L + w * n - 
    0.2e1 * L - w - 0.2e1 * xa) / (-g1 + d)) * n * w / 0.2e1;
}

double GEM_Efficiencies_v3::Fg_cath_mu2(int n, double xa, double xb)
{
    return -d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * 
    pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + 0.4e1 * 
    pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 
    0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 
    0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * 
    w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xa + 
    pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - g1 * log(pow(L, 0.2e1) * 
    pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * 
    pow(w, 0.2e1) + 0.4e1 * n * w * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * 
    xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * 
    w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 
    0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) + 
    0.4e1 * n * w * xa + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xa + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xa - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xa + pow(d, 0.2e1) + 
    0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xa + 0.4e1 * xa * xa) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 
    0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 
    0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - 
    g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w - 0.4e1 * 
    L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * 
    w * xb + 0.4e1 * xb * xb) / 0.4e1 - d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * 
    xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 
    0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + 
    0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w - 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * 
    w * xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * 
    xb + 0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 
    0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * 
    pow(w, 0.2e1) - 0.4e1 * pow(L, 0.2e1) * n - 0.6e1 * L * n * w - 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) - 0.4e1 * n * w * xb + 
    0.4e1 * pow(L, 0.2e1) + 0.4e1 * L * w + 0.8e1 * L * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) + 0.4e1 * w * 
    xb + 0.4e1 * xb * xb) / 0.4e1 + d * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + pow(n, 0.2e1) * pow(w, 0.2e1) - 
    0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) + 0.2e1 * d * g1 + pow(g1, 0.2e1) + 
    pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 + g1 * log(pow(L, 0.2e1) * pow(n, 0.2e1) + 0.2e1 * L * pow(n, 0.2e1) * w + 
    pow(n, 0.2e1) * pow(w, 0.2e1) - 0.2e1 * L * n * w + 0.4e1 * L * n * xb - 0.2e1 * n * pow(w, 0.2e1) + 0.4e1 * n * w * xb + pow(d, 0.2e1) + 
    0.2e1 * d * g1 + pow(g1, 0.2e1) + pow(w, 0.2e1) - 0.4e1 * w * xb + 0.4e1 * xb * xb) / 0.4e1 - L * atan((n * L + w * n - 0.2e1 * 
    L - w - 0.2e1 * xb) / (d + g1)) * n / 0.2e1 + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) * n / 0.2e1 - L * 
    atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * n / 0.2e1 + L * atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * n / 0.2e1 + L * 
    atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * n / 0.2e1 - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / 
    (d + g1)) * n / 0.2e1 + L * atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * n / 0.2e1 - L * atan((n * L + w * n - w + 0.2e1 * xb) / 
    (d + g1)) * n / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) * n * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * 
    L - w + 0.2e1 * xb) / (d + g1)) * n * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 + atan((n * L + w * 
    n - w + 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 + L * 
    atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) - L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) - L * 
    atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) + L * atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) + 
    atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xb) / (d + g1)) * 
    xb - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xb) / 
    (d + g1)) * xb + atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xa) / (d + g1)) * xa - 
    atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * w / 0.2e1 + atan((n * L + w * n - w + 0.2e1 * xa) / (d + g1)) * xa - atan((n * L + 
    w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w - 0.2e1 * xa) / (d + g1)) * xa + 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * 
    xa - atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * xb + 
    atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * w / 0.2e1 - atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * xb - 
    atan((n * L + w * n - 0.2e1 * L - w + 0.2e1 * xa) / (d + g1)) * n * w / 0.2e1 + atan((n * L + w * n - w - 0.2e1 * xb) / (d + g1)) * 
    n * w / 0.2e1 - atan((n * L + w * n - w + 0.2e1 * xb) / (d + g1)) * n * w / 0.2e1;
}

double GEM_Efficiencies_v3::lambda_cath_g(double xa, double xb)
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += Fg_cath_lambda(n, xa, xb);
    }
    
    result += fg_cath_lambda(xa, xb);
    
    return result;    
}

double GEM_Efficiencies_v3::mu1_cath_g(double xa, double xb)
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += Fg_cath_mu1(n, xa, xb);
    }
    
    result += fg_cath_mu1(xa, xb);
    
    return result;   
}

double GEM_Efficiencies_v3::mu2_cath_g(double xa, double xb)
{
    double result = 0.0;
    
    for ( int n = 2; n <= N; n++ )
    {
	result += Fg_cath_mu2(n, xa, xb);
    }
    
    result += fg_cath_mu2(xa, xb);
    
    return result;   
}

