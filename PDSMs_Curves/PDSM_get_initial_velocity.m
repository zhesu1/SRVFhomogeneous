function [v, qv] = PDSM_get_initial_velocity(f1,q1,f2,q2)

v=invRieExpOnSLn(f1\f2);
qv = q2 - q1;
