TITLE Mod file for component: Component(id=Im type=ionChannelHH)

COMMENT

    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)
         org.neuroml.export  v1.5.3
         org.neuroml.model   v1.5.3
         jLEMS               v0.9.9.0

ENDCOMMENT

NEURON {
    SUFFIX Im
    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!
    
    RANGE gion                           
    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!
    RANGE conductance                       : parameter
    
    RANGE g                                 : exposure
    
    RANGE fopen                             : exposure
    RANGE m_instances                       : parameter
    
    RANGE m_alpha                           : exposure
    
    RANGE m_beta                            : exposure
    
    RANGE m_tau                             : exposure
    
    RANGE m_inf                             : exposure
    
    RANGE m_rateScale                       : exposure
    
    RANGE m_fcond                           : exposure
    RANGE m_forwardRate_rate                : parameter
    RANGE m_forwardRate_midpoint            : parameter
    RANGE m_forwardRate_scale               : parameter
    
    RANGE m_forwardRate_r                   : exposure
    RANGE m_reverseRate_rate                : parameter
    RANGE m_reverseRate_midpoint            : parameter
    RANGE m_reverseRate_scale               : parameter
    
    RANGE m_reverseRate_r                   : exposure
    RANGE m_q10Settings_fixedQ10            : parameter
    
    RANGE m_q10Settings_q10                 : exposure
    RANGE conductanceScale                  : derived variable
    RANGE fopen0                            : derived variable
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!
    
    conductance = 1.0E-5 (uS)
    m_instances = 1 
    m_forwardRate_rate = 0.0033000002 (kHz)
    m_forwardRate_midpoint = -35 (mV)
    m_forwardRate_scale = 10 (mV)
    m_reverseRate_rate = 0.0033000002 (kHz)
    m_reverseRate_midpoint = -35 (mV)
    m_reverseRate_scale = -10 (mV)
    m_q10Settings_fixedQ10 = 2.9528825 
}

ASSIGNED {
    
    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel
    v (mV)
    celsius (degC)
    temperature (K)
    ek (mV)
    ik (mA/cm2)
    
    
    m_forwardRate_r (kHz)                  : derived variable
    
    m_reverseRate_r (kHz)                  : derived variable
    
    m_q10Settings_q10                      : derived variable
    
    m_rateScale                            : derived variable
    
    m_alpha (kHz)                          : derived variable
    
    m_beta (kHz)                           : derived variable
    
    m_fcond                                : derived variable
    
    m_inf                                  : derived variable
    
    m_tau (ms)                             : derived variable
    
    conductanceScale                       : derived variable
    
    fopen0                                 : derived variable
    
    fopen                                  : derived variable
    
    g (uS)                                 : derived variable
    rate_m_q (/ms)
    
}

STATE {
    m_q  
    
}

INITIAL {
    ek = -85.0
    
    temperature = celsius + 273.15
    
    rates()
    rates() ? To ensure correct initialisation.
    
    m_q = m_inf
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    
    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=Im type=ionChannelHH), from conductanceScaling; null
    ? Path not present in component, using factor: 1
    
    conductanceScale = 1 
    
    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=Im type=ionChannelHH), from gates; Component(id=m type=gateHHrates)
    ? multiply applied to all instances of fcond in: <gates> ([Component(id=m type=gateHHrates)]))
    fopen0 = m_fcond ? path based, prefix = 
    
    fopen = conductanceScale  *  fopen0 ? evaluable
    g = conductance  *  fopen ? evaluable
    gion = gmax * fopen 
    
    ik = gion * (v - ek)
    
}

DERIVATIVE states {
    rates()
    m_q' = rate_m_q 
    
}

PROCEDURE rates() {
    
    m_forwardRate_r = m_forwardRate_rate  * exp((v -  m_forwardRate_midpoint )/ m_forwardRate_scale ) ? evaluable
    m_reverseRate_r = m_reverseRate_rate  * exp((v -  m_reverseRate_midpoint )/ m_reverseRate_scale ) ? evaluable
    m_q10Settings_q10 = m_q10Settings_fixedQ10 ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=m type=gateHHrates), from q10Settings; Component(id=null type=q10Fixed)
    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10Fixed)]))
    m_rateScale = m_q10Settings_q10 ? path based, prefix = m_
    
    ? DerivedVariable is based on path: forwardRate/r, on: Component(id=m type=gateHHrates), from forwardRate; Component(id=null type=HHExpRate)
    m_alpha = m_forwardRate_r ? path based, prefix = m_
    
    ? DerivedVariable is based on path: reverseRate/r, on: Component(id=m type=gateHHrates), from reverseRate; Component(id=null type=HHExpRate)
    m_beta = m_reverseRate_r ? path based, prefix = m_
    
    m_fcond = m_q ^ m_instances ? evaluable
    m_inf = m_alpha /( m_alpha + m_beta ) ? evaluable
    m_tau = 1/(( m_alpha + m_beta ) *  m_rateScale ) ? evaluable
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    
     
    rate_m_q = ( m_inf  -  m_q ) /  m_tau ? Note units of all quantities used here need to be consistent!
    
     
    
     
    
     
    
     
    
}
