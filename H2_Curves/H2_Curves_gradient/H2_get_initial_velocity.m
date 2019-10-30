function [v, qv] = H2_get_initial_velocity(f1,q1,f2,q2)

v=H2_invRieExpOnSLn(f1\f2);
qv = q2 - q1;
