function [pars_rightward,pars_leftward,pars_NoGo,percent_right,percent_left,percent_NoGo,...
    ste_right,ste_left,ste_NoGo,trials_per_cond,min_c,max_c] = psych_curve_nogo(RM,Stimuli,unique_conditions,n_conditions)

%% function that calculates input for psychometric curve

% RespMade              should be the response made; ie going left or right
% conditions            should be the contrasts or amplitudes, in a single vector form, ie [-0.5, -0.2, 0.2, 0.5] where - denotes left/high, + right/low
%                       and that has length n_trials; so it is a list of sequences of the contrasts/amplitudes
% n_trials              should be the number of total trials, excluding baits

% unique_conditions       = unique(Stimuli);   -> should be number of unique conditions, ie how many contrasts/amplitudes were given in total
% n_conditions            = length(unique_conditions);                     % number of conditions in total, so if have 3 contrasts on each side, above vector would have length 2x3, but this is a scalar that has value 2x3

% written by Elina Jacobs, ULC Cortexlab

percent_right    = zeros(1,n_conditions);
percent_left     = zeros(1,n_conditions);
percent_NoGo     = zeros(1,n_conditions);
trials_per_cond  = zeros(1,n_conditions);

% in RespMade: -1 turn right for high freq/left contrast, 1 turn left for low freq/right contrast, 0 for NoGo
RespMade = RM;
RespMade(find(RespMade==1))=0;
RespMade = -RespMade;   %right turns are 1, everything else is zero
for ic = 1:n_conditions
    my_response = RespMade( Stimuli == unique_conditions(ic) ); % dc = diffcontrasts == contrasts_given(ic) returns a vector of length trials_per_cond(ic) where there is a 1 at the current contrast, and a 0 everywhere else, r( ) returns the turn, left or right, at those given contrast trials
    percent_right(ic) = (mean(my_response));
    trials_per_cond(ic) = length(my_response);
    ste_right(ic) = std(my_response)/sqrt(length(my_response));
end

RespMade = RM;
RespMade(find(RespMade==-1))=0;
for ic = 1:n_conditions
    my_response = RespMade( Stimuli == unique_conditions(ic) ); % dc = diffcontrasts == contrasts_given(ic) returns a vector of length trials_per_cond(ic) where there is a 1 at the current contrast, and a 0 everywhere else, r( ) returns the turn, left or right, at those given contrast trials
    percent_left(ic) = (mean(my_response));
    ste_left(ic) = std(my_response)/sqrt(length(my_response));
end

RespMade = RM;
RespMade(find(RespMade==1))=-1;
RespMade(find(RespMade==0))=1;
RespMade(find(RespMade==-1))=0; %NoGo trials become 1, everything else zero
for ic = 1:n_conditions
    my_response = RespMade( Stimuli == unique_conditions(ic) ); % dc = diffcontrasts == contrasts_given(ic) returns a vector of length trials_per_cond(ic) where there is a 1 at the current contrast, and a 0 everywhere else, r( ) returns the turn, left or right, at those given contrast trials
    percent_NoGo(ic) = (mean(my_response));
    ste_NoGo(ic) = std(my_response)/sqrt(length(my_response));
end

min_c = min(unique_conditions);
max_c = max(unique_conditions);

parstart = [ 0, 1, 0.05, 0.05 ];
parmin   = [min_c 0 0 0];
parmax   = [max_c 30 0.40 0.40];

pars_rightward = mle_fit_psycho([unique_conditions;trials_per_cond;percent_right], 'erf_psycho_2gammas', parstart,parmin,parmax);
pars_leftward = mle_fit_psycho([unique_conditions;trials_per_cond;percent_left], 'erf_psycho_2gammas', parstart,parmin,parmax);
pars_NoGo = mle_fit_psycho([unique_conditions;trials_per_cond;percent_NoGo], 'erf_psycho_2gammas', parstart,parmin,parmax);

