/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__SK_E2
#define _nrn_initial _nrn_initial__SK_E2
#define nrn_cur _nrn_cur__SK_E2
#define _nrn_current _nrn_current__SK_E2
#define nrn_jacob _nrn_jacob__SK_E2
#define nrn_state _nrn_state__SK_E2
#define _net_receive _net_receive__SK_E2 
#define rates rates__SK_E2 
#define states states__SK_E2 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define conductance _p[1]
#define z_instances _p[2]
#define z_timeCourse_TIME_SCALE _p[3]
#define z_timeCourse_VOLT_SCALE _p[4]
#define z_timeCourse_CONC_SCALE _p[5]
#define z_steadyState_TIME_SCALE _p[6]
#define z_steadyState_VOLT_SCALE _p[7]
#define z_steadyState_CONC_SCALE _p[8]
#define gion _p[9]
#define z_timeCourse_V _p[10]
#define z_timeCourse_ca_conc _p[11]
#define z_timeCourse_t _p[12]
#define z_steadyState_V _p[13]
#define z_steadyState_ca_conc _p[14]
#define z_steadyState_x _p[15]
#define z_rateScale _p[16]
#define z_fcond _p[17]
#define z_inf _p[18]
#define z_tauUnscaled _p[19]
#define z_tau _p[20]
#define conductanceScale _p[21]
#define fopen0 _p[22]
#define fopen _p[23]
#define g _p[24]
#define z_q _p[25]
#define temperature _p[26]
#define ek _p[27]
#define ik _p[28]
#define cai _p[29]
#define cao _p[30]
#define rate_z_q _p[31]
#define Dz_q _p[32]
#define v _p[33]
#define _g _p[34]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ik	*_ppvar[2]._pval
#define _ion_dikdv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_SK_E2", _hoc_setdata,
 "rates_SK_E2", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_SK_E2", "S/cm2",
 "conductance_SK_E2", "uS",
 "z_timeCourse_TIME_SCALE_SK_E2", "ms",
 "z_timeCourse_VOLT_SCALE_SK_E2", "mV",
 "z_timeCourse_CONC_SCALE_SK_E2", "mM",
 "z_steadyState_TIME_SCALE_SK_E2", "ms",
 "z_steadyState_VOLT_SCALE_SK_E2", "mV",
 "z_steadyState_CONC_SCALE_SK_E2", "mM",
 "gion_SK_E2", "S/cm2",
 "z_timeCourse_t_SK_E2", "ms",
 "z_tauUnscaled_SK_E2", "ms",
 "z_tau_SK_E2", "ms",
 "g_SK_E2", "uS",
 0,0
};
 static double delta_t = 0.01;
 static double z_q0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"SK_E2",
 "gmax_SK_E2",
 "conductance_SK_E2",
 "z_instances_SK_E2",
 "z_timeCourse_TIME_SCALE_SK_E2",
 "z_timeCourse_VOLT_SCALE_SK_E2",
 "z_timeCourse_CONC_SCALE_SK_E2",
 "z_steadyState_TIME_SCALE_SK_E2",
 "z_steadyState_VOLT_SCALE_SK_E2",
 "z_steadyState_CONC_SCALE_SK_E2",
 0,
 "gion_SK_E2",
 "z_timeCourse_V_SK_E2",
 "z_timeCourse_ca_conc_SK_E2",
 "z_timeCourse_t_SK_E2",
 "z_steadyState_V_SK_E2",
 "z_steadyState_ca_conc_SK_E2",
 "z_steadyState_x_SK_E2",
 "z_rateScale_SK_E2",
 "z_fcond_SK_E2",
 "z_inf_SK_E2",
 "z_tauUnscaled_SK_E2",
 "z_tau_SK_E2",
 "conductanceScale_SK_E2",
 "fopen0_SK_E2",
 "fopen_SK_E2",
 "g_SK_E2",
 0,
 "z_q_SK_E2",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 35, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	conductance = 1e-05;
 	z_instances = 1;
 	z_timeCourse_TIME_SCALE = 1;
 	z_timeCourse_VOLT_SCALE = 1;
 	z_timeCourse_CONC_SCALE = 1e+06;
 	z_steadyState_TIME_SCALE = 1;
 	z_steadyState_VOLT_SCALE = 1;
 	z_steadyState_CONC_SCALE = 1e+06;
 	_prop->param = _p;
 	_prop->param_size = 35;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 prop_ion = need_memb(_k_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _SK_E2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", 2.0);
 	ion_reg("k", 1.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 35, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SK_E2 /home/vrhaynes/workspace/Haynes2021_EAPs/test_bbp_models/NMLCL000637-NEURON/x86_64/SK_E2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Mod file for component: Component(id=SK_E2 type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dz_q = rate_z_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dz_q = Dz_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    z_q = z_q - dt*(- ( rate_z_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
   double _lcaConc ;
 _lcaConc = cai ;
   z_timeCourse_V = v / z_timeCourse_VOLT_SCALE ;
   z_timeCourse_ca_conc = _lcaConc / z_timeCourse_CONC_SCALE ;
   z_timeCourse_t = 1.0 * z_timeCourse_TIME_SCALE ;
   z_steadyState_V = v / z_steadyState_VOLT_SCALE ;
   z_steadyState_ca_conc = _lcaConc / z_steadyState_CONC_SCALE ;
   z_steadyState_x = 1.0 / ( 1.0 + pow( ( 4.3e-10 / z_steadyState_ca_conc ) , 4.8 ) ) ;
   z_rateScale = 1.0 ;
   z_fcond = pow( z_q , z_instances ) ;
   z_inf = z_steadyState_x ;
   z_tauUnscaled = z_timeCourse_t ;
   z_tau = z_tauUnscaled / z_rateScale ;
   rate_z_q = ( z_inf - z_q ) / z_tau ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  z_q = z_q0;
 {
   ek = - 85.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   z_q = z_inf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   conductanceScale = 1.0 ;
   fopen0 = z_fcond ;
   fopen = conductanceScale * fopen0 ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ik = gion * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(z_q) - _p;  _dlist1[0] = &(Dz_q) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/vrhaynes/workspace/Haynes2021_EAPs/test_bbp_models/NMLCL000637-NEURON/SK_E2.mod";
static const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=SK_E2 type=ionChannelHH)\n"
  "\n"
  "COMMENT\n"
  "\n"
  "    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)\n"
  "         org.neuroml.export  v1.5.3\n"
  "         org.neuroml.model   v1.5.3\n"
  "         jLEMS               v0.9.9.0\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX SK_E2\n"
  "    USEION ca READ cai,cao VALENCE 2\n"
  "    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE z_instances                       : parameter\n"
  "    \n"
  "    RANGE z_tau                             : exposure\n"
  "    \n"
  "    RANGE z_inf                             : exposure\n"
  "    \n"
  "    RANGE z_rateScale                       : exposure\n"
  "    \n"
  "    RANGE z_fcond                           : exposure\n"
  "    RANGE z_timeCourse_TIME_SCALE           : parameter\n"
  "    RANGE z_timeCourse_VOLT_SCALE           : parameter\n"
  "    RANGE z_timeCourse_CONC_SCALE           : parameter\n"
  "    \n"
  "    RANGE z_timeCourse_t                    : exposure\n"
  "    RANGE z_steadyState_TIME_SCALE          : parameter\n"
  "    RANGE z_steadyState_VOLT_SCALE          : parameter\n"
  "    RANGE z_steadyState_CONC_SCALE          : parameter\n"
  "    \n"
  "    RANGE z_steadyState_x                   : exposure\n"
  "    RANGE z_timeCourse_V                    : derived variable\n"
  "    RANGE z_timeCourse_ca_conc              : derived variable\n"
  "    RANGE z_steadyState_V                   : derived variable\n"
  "    RANGE z_steadyState_ca_conc             : derived variable\n"
  "    RANGE z_tauUnscaled                     : derived variable\n"
  "    RANGE conductanceScale                  : derived variable\n"
  "    RANGE fopen0                            : derived variable\n"
  "    \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    \n"
  "    (nA) = (nanoamp)\n"
  "    (uA) = (microamp)\n"
  "    (mA) = (milliamp)\n"
  "    (A) = (amp)\n"
  "    (mV) = (millivolt)\n"
  "    (mS) = (millisiemens)\n"
  "    (uS) = (microsiemens)\n"
  "    (molar) = (1/liter)\n"
  "    (kHz) = (kilohertz)\n"
  "    (mM) = (millimolar)\n"
  "    (um) = (micrometer)\n"
  "    (umol) = (micromole)\n"
  "    (S) = (siemens)\n"
  "    \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    \n"
  "    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!\n"
  "    \n"
  "    conductance = 1.0E-5 (uS)\n"
  "    z_instances = 1 \n"
  "    z_timeCourse_TIME_SCALE = 1 (ms)\n"
  "    z_timeCourse_VOLT_SCALE = 1 (mV)\n"
  "    z_timeCourse_CONC_SCALE = 1000000 (mM)\n"
  "    z_steadyState_TIME_SCALE = 1 (ms)\n"
  "    z_steadyState_VOLT_SCALE = 1 (mV)\n"
  "    z_steadyState_CONC_SCALE = 1000000 (mM)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    ek (mV)\n"
  "    ik (mA/cm2)\n"
  "    \n"
  "    cai (mM)\n"
  "    \n"
  "    cao (mM)\n"
  "    \n"
  "    \n"
  "    z_timeCourse_V                         : derived variable\n"
  "    \n"
  "    z_timeCourse_ca_conc                   : derived variable\n"
  "    \n"
  "    z_timeCourse_t (ms)                    : derived variable\n"
  "    \n"
  "    z_steadyState_V                        : derived variable\n"
  "    \n"
  "    z_steadyState_ca_conc                  : derived variable\n"
  "    \n"
  "    z_steadyState_x                        : derived variable\n"
  "    \n"
  "    z_rateScale                            : derived variable\n"
  "    \n"
  "    z_fcond                                : derived variable\n"
  "    \n"
  "    z_inf                                  : derived variable\n"
  "    \n"
  "    z_tauUnscaled (ms)                     : derived variable\n"
  "    \n"
  "    z_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopen0                                 : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_z_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    z_q  \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    ek = -85.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    z_q = z_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=SK_E2 type=ionChannelHH), from conductanceScaling; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    conductanceScale = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=SK_E2 type=ionChannelHH), from gates; Component(id=z type=gateHHtauInf)\n"
  "    ? multiply applied to all instances of fcond in: <gates> ([Component(id=z type=gateHHtauInf)]))\n"
  "    fopen0 = z_fcond ? path based, prefix = \n"
  "    \n"
  "    fopen = conductanceScale  *  fopen0 ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ik = gion * (v - ek)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    z_q' = rate_z_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    LOCAL caConc\n"
  "    \n"
  "    caConc = cai\n"
  "    \n"
  "    z_timeCourse_V = v /  z_timeCourse_VOLT_SCALE ? evaluable\n"
  "    z_timeCourse_ca_conc = caConc /  z_timeCourse_CONC_SCALE ? evaluable\n"
  "    z_timeCourse_t = 1.0  *  z_timeCourse_TIME_SCALE ? evaluable\n"
  "    z_steadyState_V = v /  z_steadyState_VOLT_SCALE ? evaluable\n"
  "    z_steadyState_ca_conc = caConc /  z_steadyState_CONC_SCALE ? evaluable\n"
  "    z_steadyState_x = 1/(1+(4.3e-10/ z_steadyState_ca_conc )^4.8) ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=z type=gateHHtauInf), from q10Settings; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    z_rateScale = 1 \n"
  "    \n"
  "    z_fcond = z_q ^ z_instances ? evaluable\n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=z type=gateHHtauInf), from steadyState; Component(id=null type=SK_E2_z_inf_inf)\n"
  "    z_inf = z_steadyState_x ? path based, prefix = z_\n"
  "    \n"
  "    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=z type=gateHHtauInf), from timeCourse; Component(id=null type=SK_E2_z_tau_tau)\n"
  "    z_tauUnscaled = z_timeCourse_t ? path based, prefix = z_\n"
  "    \n"
  "    z_tau = z_tauUnscaled  /  z_rateScale ? evaluable\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_z_q = ( z_inf  -  z_q ) /  z_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "}\n"
  "\n"
  ;
#endif
