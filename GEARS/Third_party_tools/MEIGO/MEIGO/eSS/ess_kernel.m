
function [Results]=ess_kernel(problem,opts,varargin)
% function [Results]=ess_kernel(problem,opts,varargin)
% Function   : eSS Release 2016A - last revised Nov-2016 - JEega, JRB
% Written by : Process Engineering Group IIM-CSIC (jegea@iim.csic.es)
% Created on : 12/03/2008
% Last Update: 11/11/2016 
% Email      : gingproc@iim.csic.es
%
% (c) CSIC, Spanish Council for Scientific Research
%
% Global optimization algorithm for MINLP's based on Scatter Search
%
%   eSS attempts to solve problems of the form:
%       min F(x)  subject to:  ceq(x) = 0 (equality constraints)
%        x                     c_L <= c(x) <= c_U (inequality constraints)
%                              x_L <= x <= x_U (bounds on the decision
%                              variables)
%
% Constraint functions, if applicable, must be declared in the same script as the
% objective function as a second output argument:
% e.g., [f,c]=myfunction(x)
%       f=objective function value
%       g=vector containing the constraints
%       return
% 
% DW: [f, c, gf, gc] = myfunction(x)
%       gf, gc: gradient of objective function and constraints
%
% Please have a look at the manual before using eSS
%
%       USAGE: Results=ess_kernel(problem,opts,p1,p2,....,pn);
%
% INPUT PARAMETERS:
% ****************
%   problem - Structure containing problem settings
%       problem.f   = Name of the file containing the objective
%                     function (String)
%       problem.x_L = Lower bounds of decision variables (vector)
%       problem.x_U = Upper bounds of decision variables (vector)
%       problem.x_0 = Initial point(s) (optional; vector or matrix)
%       problem.f_0 = Function values of initial point(s) (optional) These
%       value MUST correspond to feasible points.
%
%       NOTE:The dimension of f_0 and x_0 may be different. For example, if
%       we want to introduce 5 initial points but we only know the values
%       for 3 of them, x_0 would have 5 rows whereas f_0 would have only 3
%       elements. It is mandatory that the first 3 rows of x_0 correspond
%       to the values of f_0
%
%       Additionally, fill the following fields if your problem has
%       non-linear constraints
%       problem.neq     = Number of equality constraints (Integer; do not define it
%                         if there are no equality constraints)
%       problem.c_L     = Lower bounds of nonlinear inequality constraints
%                        (vector)
%       problem.c_U     = Upper bounds of nonlinear inequality constraints
%                         (vector)
%       problem.int_var = Number of integer variables (Integer)
%       problem.bin_var = Number of binary variables  (Integer)
%       problem.vtr     = Objective function value to be reached (optional)
%
%
% NOTE: The order of decision variables is x=[cont int bin]
%
%   opts - Structure containing options (if set as opts=[] defaults options
%          will be loaded) Type "ssm_kernel" or "ssm_kernel('defaults')" to
%          get the default options
%
%       User options
%           opts.maxeval    = Maximum number of function evaluations
%                             (Default 1000)
%           opts.maxtime    = Maximum CPU time in seconds (Default 60)
%           opts.iterprint  = Print each iteration on screen: 0-Deactivated
%                             1-Activated (Default 1)
%           opts.plot       = Plots convergence curves: 0-Deactivated,
%                             1-Plot curves on line, 2-Plot final results
%                             (Default 0)
%           opts.weight     = Weight that multiplies the penalty term added
%                             to the objective function in constrained
%                             problems (Default 1e6)
%           opts.log_var    = Indexes of the variables which will be used
%                             to generate diverse solutions in different
%                             orders of magnitude (vector)
%           opts.log_level  = Definition of what actions should be carried
%                             out in log space: 1 - only the refset, 2 - 
%                             all actions except the local search and 3 - 
%                             everything. (scaler) (Defualt 1)
%           opts.tolc       = Maximum absolute violation of the constraints
%                             (Default 1e-5)
%           opts.prob_bound = Probability (0-1) of biasing the search towards
%                             the bounds (Default 0.5)
%           opts.inter_save = Saves results in a mat file in intermediate
%                             iterations. Useful for very long runs
%                             (Binary; Default = 0)
%           opts.n_stuck    = Number of consecutive iterations without
%                             significant improvement before the search
%                             stops. Default = 0 (not applicable);
%
%       Global options
%           opts.dim_refset         = Number of elements in Refset
%                                     (Integer; automatically calculated)
%           opts.ndiverse           = Number of solutions generated by the
%                                     diversificator (Default 10*nvar)
%           opts.combination        = Type of combination of Refset
%                                     elements (Default 1)
%                                     1: hyper-rectangles
%                                     2: linear combinations
%
%       Local options
%           opts.local.solver               = Choose local solver
%                                             0: Local search deactivated
%                                             'fmincon'(Default),
%                                             'fminsearch', 'solnp'
%                                             'n2fb','dn2fb','dhc','hooke'
%                                             'ipopt','misqp','lsqnonlin'
%           opts.local.tol                  = Level of tolerance in local
%                                             search. 1: Relaxed; 2: Medium
%                                             3: Tight (default: 2 in
%                                             intermediate searches and 3
%                                             in the final stage).
%           opts.local.iterprint            = Print each iteration of local
%                                             solver on screen (Binary;
%                                             default = 0)
%           opts.local.n1                   = Number of iterations
%                                             before applying
%                                             local search for the 1st time
%                                             (Default 1)
%           opts.local.n2                   = Minimum number of iterations
%                                             in the global
%                                             phase between 2 local calls
%                                             (Default 10)
%           opts.local.balance              = Balances between quality (=0)
%                                             and diversity (=1) for
%                                             choosing initial points for
%                                             the local search (default
%                                             0.5)
%           opts.local.finish               = Applies local search to the
%                                             best solution found once the
%                                             optimization if finished
%                                             (same values as
%                                             opts.local.solver)
%           opts.local.bestx                = When activated (i.e. =1) only
%                                             applies local search to the
%                                             best solution found to
%                                             date,ignoring filters
%                                             (Default=0)
%
%
%   p1,p2... :  optional input parameters to be passed to the objective
%   function
%
%
% OUTPUT PARAMETERS:
% *****************
% A file called "ssm_report.mat" is generated containing.
%
%   problem - Structure containing problem settings
%   opts    - Structure containing all options
%   Results - Structure containing results
%
% Fields in Results
%       Results.fbest                   = Best objective function value
%                                         found after the optimization
%       Results.xbest                   = Vector providing the best
%                                         function value
%       Results.cpu_time                = Time in seconds consumed in the
%                                         optimization
%       Results.f                       = Vector containing the best
%                                         objective function value after each
%                                         iteration
%       Results.x                       = Matrix containing the best vector
%                                         after each iteration
%       Results.time                    = Vector containing the cpu time
%                                         consumed after each iteration
%       Results.neval                   = Vector containing the number of
%                                         function evaluations after each
%                                         iteration
%       Results.numeval                 = Number of function evaluations
%       Results.local_solutions         = Local solutions found by the
%                                         local solver (in rows)
%       Results.local_solutions_values  = Function values of the local
%                                         solutions
%       Results.end_crit                = Criterion to finish the
%                                         optimization
%                                         1: Maximal number of function
%                                         evaluations achieved
%                                         2: Maximum allowed CPU Time
%                                         achieved
%                                         3: Value to reach achieved
%
%              NOTE: To plot convergence curves type:
%              stairs(Results.time,Results.f) or stairs(Results.neval,Results.f)
% REFERENCES
%
% If you use eSS and publish the results, please cite the following papers:
%
% Egea, J.A., Balsa-Canto, E., Garc�a, M.S.G., Banga, J.R. (2009) Dynamic
% optimization of nonlinear processes with an enhanced scatter search
% method. Industrial & Engineering Chemistry Research 49(9): 4388-4401.
%
% Egea, J.A., Mart�, R., Banga, J.R. (2010) An evolutionary method for
% complex-process optimization. Computers & Operations Research 37(2):
% 315-324.


global input_par
global stopOptimization
stopOptimization=0;

%Initialize time
cpu_time=cputime;


% if isfield(opts,'iterprint') & opts.iterprint
%     vers='eSS R2014B  IIM-CSIC, Vigo, Spain';% AMIGO2014bench VERSION - JRB
%     fprintf('%-9s %s\n','Version:',  vers);
% end

fprintf('\n');
fprintf('------------------------------------------------------------------------------ \n');
fprintf(' eSS R2016A - Enhanced Scatter Search    \n');
fprintf('<c> IIM-CSIC, Vigo, Spain -  email: gingproc@iim.csic.es \n');
fprintf('------------------------------------------------------------------------------ \n\n');


%Initialize the number of function evaluations and a stop variable
nfuneval=0;
fin=0;
nreset=0;

if not(isfield(opts,'maxeval')) && not(isfield(opts,'maxtime'))
    fprintf('WARNING:Either opts.maxeval or opts.maxtime must be defined as a stop criterion \n')
    fprintf('Define any of these options and rerun \n')
    Results=[];
    return
else
    if not(isfield(opts,'maxeval'))
        maxeval=1e12;
    else
        maxeval=opts.maxeval;
    end
    if not(isfield(opts,'maxtime'))
        maxtime=1e12;
    else
        maxtime=opts.maxtime;
    end
end
   
%Load default values for the options
default=ssm_defaults;
loc=default.local;

%Display default options
if (nargin < 1) | (isstr(problem) & isequal(problem, 'defaults')) % pass default options
    disp('Default options returned (type "help ssm_kernel" for help).');
    fprintf('\n');
    disp('Global options')
    disp(default)
    disp('Local options')
    disp(loc)
    Results = 'End of default options display';
    return
end

if nargin<2, opts=[]; end

% check regularization for iterative display of the reg. cost and penalty
reg_flag = 0;
if nargin == 5
    if isfield(varargin{1},'nlpsol') && isfield(varargin{1}.nlpsol,'regularization') && varargin{1}.nlpsol.regularization.ison
        reg_flag = 1;
        amigo_inputs = varargin{1};
        amigo_privstruct = varargin{2};
    end
end

opts=ssm_optset(default,opts);

log_level = opts.log_level;

%Set all options

%Transform the objective function into a handle
if isa(problem.f,'function_handle')
    fobj = problem.f;
    problem.f = func2str(problem.f);
else
    fobj=str2func(problem.f);
end

isFjacDefined = 0;
if isfield(problem,'fjac') && ~isempty(problem.fjac)
    isFjacDefined = 1;
    
	if ischar(problem.fjac)
		fjac = str2func(problem.fjac);
	elseif isa(problem.fjac,'function_handle')
		fjac = problem.fjac;
	else
		error('fjac is neither a string nor a function handle.')
	end
end

log_var=opts.log_var;

%Extra input parameters for fsqp, n2fb and nomad
if nargin>2
    if strcmp(opts.local.finish,'n2fb') | strcmp(opts.local.finish,'dn2fb') |...
            strcmp(opts.local.solver,'dn2fb') | strcmp(opts.local.solver,'n2fb') |...
            strcmp(opts.local.solver,'fsqp') |strcmp(opts.local.finish,'fsqp')|...
            strcmp(opts.local.solver,'nomad') |strcmp(opts.local.finish,'nomad')|...
            strcmp(opts.local.solver,'nl2sol')|strcmp(opts.local.finish,'nl2sol')|...
            strcmp(opts.local.solver,'lbfgsb')|strcmp(opts.local.finish,'lbfgsb')
            %%% Jake set up passing variables into the cost interfaces for log level 3 % Start
            input_par = varargin;
            
                if log_level == 3
                    
                input_par{end+1} = log_var;    
                
                input_par{end+1} = fobj;
                
                    if isFjacDefined
                        
                    input_par{end+1} = fjac;
                    
                    end
                    
                end
                   
            %%% Jake set up passing variables into the cost interfaces for log level 3 % End
    end
end

x_U=problem.x_U;
x_L=problem.x_L;

%Check if bounds have the same dimension
if length(x_U)~=length(x_L)
    disp('Upper and lower bounds have different dimension!!!!')
    disp('EXITING')
    Results=[];
    return
else
    %Number of decision variables
    nvar=length(x_L);
end

if not(isfield(problem,'x_0')), x_0=[]; else, x_0=problem.x_0; end

if not(isfield(problem,'f_0'))
    f_0=[];
    fbest=inf;
else
    f_0=problem.f_0;
    %Transform in column vector
    if size(f_0,1)==1
        f_0=f_0';
    end
    [fbest,iii]=min(f_0);
    xbest=x_0(iii,:);
end
if not(isfield(problem,'neq')) | isempty(problem.neq), neq=0; else, neq=problem.neq; end
if not(isfield(problem,'c_U')), c_U=[]; c_L=[]; else, c_U=problem.c_U; c_L=problem.c_L; end
if not(isfield(problem,'int_var')) | isempty(problem.int_var), int_var=0; else, int_var=problem.int_var; end
if not(isfield(problem,'bin_var')) | isempty(problem.bin_var), bin_var=0; else, bin_var=problem.bin_var; end
if not(isfield(problem,'vtr')), vtr=[]; else, vtr=problem.vtr; end


%Change to rows if necessary
if size(x_L,2)==1, x_L=x_L'; end
if size(x_U,2)==1, x_U=x_U'; end
if size(c_L,2)==1, c_L=c_L'; end
if size(c_U,2)==1, c_U=c_U'; end

%User options
iterprint=opts.iterprint;
prob_bound=opts.prob_bound;
inter_save=opts.inter_save;

plot_results=opts.plot;
weight=opts.weight;
tolc=opts.tolc;

%Global options
dim_refset=opts.dim_refset;
if ischar(dim_refset) || isempty(dim_refset)
    nnn=roots([1  -1 -10*nvar/1]);
    iii=find(nnn>0);
    dim_refset=ceil(nnn(iii));
    if mod(dim_refset,2), dim_refset=dim_refset+1; end
    if iterprint
        fprintf('%-9s %i\n','Refset size automatically calculated:',  dim_refset);
    end
end

ndiverse=opts.ndiverse;

if ischar(ndiverse)
    
    ndiverse=10*nvar;
    if iterprint
        fprintf('%-9s %i\n','Number of diverse solutions automatically calculated:',  ndiverse);
    end
    
elseif(isempty(ndiverse))
    
    ndiverse=10*nvar;
    if iterprint
        fprintf('%-9s %i\n','Number of diverse solutions automatically calculated:',  ndiverse);
    end
end

if ndiverse<dim_refset
    ndiverse=dim_refset;
end

initiate=opts.initiate;
combin=opts.combination;
n_stuck=opts.n_stuck;

%Local options
local_solver=opts.local.solver;
local_tol=opts.local.tol;
local_iterprint=opts.local.iterprint;
local_n1=opts.local.n1;
local_n2=opts.local.n2;
if isstr(local_n1), local_n1=1; end
if isstr(local_n2), local_n2=10;end


local_balance=opts.local.balance;

if isempty(opts.local.finish)
    opts.local.finish=local_solver;
    local_finish=local_solver;
else
    local_finish=opts.local.finish;
end
local_bestx=opts.local.bestx;


%If there are equality constraints
if neq
    %Set the new bounds for all the constraints
    c_L=[zeros(1,neq) c_L];
    c_U=[zeros(1,neq) c_U];
end

if size(x_0,2)~=nvar, x_0=x_0'; end

nconst=length(c_U);

%Pasamos a logaritmo
xl_log=x_L;
xu_log=x_U;

aaa= (xl_log(log_var)==0);

xl_log(log_var(aaa))=eps;

xl_log(log_var)=log(xl_log(log_var));
xu_log(log_var)=log(xu_log(log_var));

    %%% Jake convert all eSS operations into log scale (if == 2 local solver will be transformed back). % Start
    if log_level > 1
        
    x_L(log_var(aaa))=eps;
        
    x_L(log_var)=log(x_L(log_var));

    x_U(log_var)=log(x_U(log_var));    
                
    end
    %%% Jake convert all eSS operations into log scale (if == 2 local solver will be transformed back). % End

%Check if the objecttive function has 2 output arguments in constrained
%problems
if nconst
    n_out_f=nargout(problem.f);
    if n_out_f<2
        fprintf('For constrained problems the objective function must have at least 2 output arguments \n');
        fprintf('EXITING \n');
        Results=[];
        return
    end
    if local_solver & strcmp(local_solver,'fminsearch') | strcmp(local_solver,'hooke') | strcmp(local_solver,'dhc')
        fprintf('Warning! %s might not provide the optimal solution for constrained problems \n',local_solver);
    end
end

if (int_var+bin_var) & local_solver & not(strcmp(local_solver,'misqp'))
    fprintf('For problems with integer and/or binary variables you must use MISQP as a local solver \n');
    fprintf('EXITING \n');
    Results=[];
    return
end

if local_solver & strcmp(local_solver,'n2fb') | strcmp(local_solver,'dn2fb')| strcmp(local_solver,'nl2sol')
    n_out_f=nargout(problem.f);
    if n_out_f<3
        fprintf('%s requires 3 output arguments \n',local_solver);
        fprintf('EXITING \n');
        Results=[];
        return
    else
        %Generate a random point within the bounds
        randx=rand(1,nvar).*(x_U-x_L)+x_L;
            %%% Jake evaluate the value and not log of the value % Start
            if log_level > 1
                
            randx(log_var) = exp(randx(log_var));
            
            [f g R]=feval(fobj,randx,varargin{:});

            randx(log_var) = log(randx(log_var));    
                        
            else
                
            [f g R]=feval(fobj,randx,varargin{:});
   
            end
            %%% Jake evaluate the value and not log of the value % End
        ndata=length(R);
        
        if ~isempty(g) & g
            fprintf('%s can not solve constrained problems \n',local_solver);
            fprintf('EXITING \n');
            Results=[];
            return
        end
        
        if length(R)<=1
            fprintf('The third output argument must be a vector, not a scalar \n');
            fprintf('EXITING \n');
            Results=[];
            return
        end
        nfuneval=nfuneval+1;
    end
else
    ndata=[];
end

%DW 
if opts.local.use_gradient_for_finish
    if ~local_finish | ~strcmp(local_finish,'fmincon')
        fprintf(['opts.local.use_gradient_for_finish == 1, can only be used with fmincon for the moment\n']);
        fprintf('EXITING \n');
        Results=[];
    end
    if nargout(problem.f) < 3
        fprintf(['opts.local.use_gradient_for_finish == 1, but no gradients are provided by "' problem.f '"\n']);
        fprintf('EXITING \n');
        Results=[];
        return
    end
end

if local_finish & isempty(ndata) & strcmp(local_finish,'n2fb') | strcmp(local_finish,'dn2fb') | strcmp(local_finish,'nl2sol')
    n_out_f=nargout(problem.f);
    if n_out_f<3
        fprintf('%s requires 3 output arguments \n',local_solver);
        fprintf('EXITING \n');
        Results=[];
        return
    else
        %Generate a random point within the bounds
        randx=rand(1,nvar).*(x_U-x_L)+x_L;
        
            if log_level > 1
                
            randx(log_var) = exp(randx(log_var));
            
            [f g R]=feval(fobj,randx,varargin{:});

            randx(log_var) = log(randx(log_var));    
                        
            else
                
            [f g R]=feval(fobj,randx,varargin{:});
   
            end
        ndata=length(R);
        if ~isempty(g) & g
            fprintf('%s can not solve constrained problems \n',local_solver);
            fprintf('EXITING \n');
            Results=[];
            return
        end
        
        if length(R)<=1
            fprintf('The third output argument must be a vector, not a scalar \n');
            fprintf('EXITING \n');
            Results=[];
            return
        end
        nfuneval=nfuneval+1;
    end
end

%Build auxiliary files in case local search is activated
if local_solver
    if isFjacDefined
        ssm_aux_local(problem.f,problem.fjac,x_L,x_U,c_L,c_U,neq,local_solver,nvar, log_level, log_var, varargin{:});
    else
        ssm_aux_local(problem.f,[],x_L,x_U,c_L,c_U,neq,local_solver,nvar, log_level, log_var, varargin{:});
    end
end
% ssm_aux_local(problem.f,x_L,x_U,c_L,c_U,neq,local_solver,nvar,varargin{:});
%Initialize variables
stage_1=1;
n_minimo=0;
n_critico=local_n1;
denom=2;
suit_to_local=ones(1,dim_refset);

local_solutions=[];
local_solutions_values=[];
initial_points=[];

Results.f=[];
Results.x=[];
Results.time=[];
Results.neval=[];
xbest=[];

%Possible combinations among the dim_refset elements in Refset, taken in pairs
ncomb=nchoosek(1:dim_refset,2);
MaxSubSet=(dim_refset^2-dim_refset)/2;
MaxSubSet2=2*MaxSubSet;


%We generate 5 solutions
solutions=zeros(ndiverse+5,nvar);

for i=1:5
    solutions(i,:)=(rand(1,nvar)+i-1)/5;
end

solutions(6:ndiverse+5,:)=rand(ndiverse,nvar);
solutions=solutions.*repmat(xu_log-xl_log,ndiverse+5,1)+repmat(xl_log,ndiverse+5,1);

% BUG FIXED - JRB - JUNE 2014
%if (int_var) | (bin_var)
%    solutions=ssm_round_int(solutions,int_var+bin_var,x_L,x_U);
%end

	if log_level < 2
	
	for i=1:ndiverse+5
		   solutions(i,log_var)=exp(solutions(i,log_var));
	end

	end

l_f_0=length(f_0);
l_x_0=size(x_0,1);

%Put x_0 without their corresponding f_0 values at the beginning of solutions
solutions=[x_0(l_f_0+1:end,:); solutions];

temp=l_x_0-l_f_0;

sol=zeros(ndiverse+5+temp,nvar);
sol_val=zeros(ndiverse+5+temp,1);
sol_val_pen=zeros(ndiverse+5+temp,1);
sol_pen=zeros(ndiverse+5+temp,1);
sol_nlc=zeros(ndiverse+5+temp,nconst);

counter=1;
%Evaluate the set of solutions

for i=1:ndiverse+5+temp
    [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(solutions(i,:),x_L,x_U,fobj,nconst,c_L,c_U,...
        tolc,weight,int_var,bin_var,log_level,log_var,varargin{:});
    nfuneval=nfuneval+1;
    if include
        sol(counter,:)=x;
        sol_val(counter)=val;
        sol_val_pen(counter)=val_penalty;
        sol_pen(counter)=pena;
        sol_nlc(counter,:)=nlc;
        counter=counter+1;
    end
    
end

sol(counter:end,:)=[];
sol_val(counter:end)=[];
sol_val_pen(counter:end)=[];
sol_pen(counter:end)=[];
sol_nlc(counter:end,:)=[];


%Add points with f_0 values. We assume feasiblity

sol=[x_0(1:l_f_0,:); sol];
sol_val=[f_0; sol_val];
sol_val_pen=[f_0; sol_val_pen];   %We assume feasible points in f_0
sol_pen=[zeros(l_f_0,1);sol_pen];         %We assume feasible points in f_0
sol_nlc=[NaN*ones(l_f_0,nconst); sol_nlc];     %We do not know these values, but we do not care. We assume feasibility

%% Update: 11/11/2016 
%% The following added by JEgea on Nov-2016 to handle REfSet with infeasible solutions
%% JRB reported the problem found in constrained problems: when initial point is infeasible, and
%% all the elements in the generated RefSet are also infeasible, eSS was printing obj=Inf, which was misleading
if isfield(problem,'x_0')
    %This is if only one solution in x_0 is defined
    if sol_pen(1) > 0
     fprintf('Initial Point (x_0): INFEASIBLE SOLUTION. Penalized value: %g \n', sol_val_pen(1));
    else
       fprintf('Initial Point (x_0): FEASIBLE SOLUTION. Objective function value: %g \n', sol_val_pen(1));
    end
 end
%%


%Initialize Refset
Refset=zeros(dim_refset,nvar);
Refset_values=zeros(dim_refset,1);
Refset_values_penalty=zeros(dim_refset,1);
Refset_nlc=[];
penalty=[];


[aaa, iii]=sort(sol_val_pen);
first_members=ceil(dim_refset/2);
last_members=dim_refset-first_members;


%Create the initial sorted Refset

Refset=sol(iii(1:first_members),:);
Refset_values=sol_val(iii(1:first_members));
Refset_values_penalty=sol_val_pen(iii(1:first_members));
penalty=sol_pen(iii(1:first_members));
Refset_nlc=sol_nlc(iii(1:first_members),:);

%The rest of Refset members are chosen randomly
iii(1:first_members)=[];
ooo=randperm(length(iii));
%ooo=randperm(ndiverse+2+size(x_0,1)-first_members);

Refset=[Refset; sol(iii(ooo(1:last_members)),:)];
Refset_values=[Refset_values; sol_val(iii(ooo(1:last_members)))];
Refset_values_penalty=[Refset_values_penalty; sol_val_pen(iii(ooo(1:last_members)))];
penalty=[penalty; sol_val_pen(iii(ooo(1:last_members)))];
Refset_nlc=[Refset_nlc; sol_nlc(iii(ooo(1:last_members)),:)];

%Check best value
ggg=find(penalty<=tolc);
if not(isempty(ggg))
    [fbest,iii]=min(Refset_values_penalty(ggg));
    xbest=Refset((ggg(iii)),:);
end


iter=0;
if iterprint==1
    if reg_flag
        amigo_privstruct.theta = xbest;
        amigo_privstruct.theta = AMIGO_scale_theta_back(amigo_privstruct.theta, amigo_inputs);
        [objReg] = AMIGO_PEcostreg(amigo_privstruct.theta,amigo_inputs,amigo_privstruct);
        reg_alpha = amigo_inputs.nlpsol.regularization.alpha;
        fprintf('Initial Pop: NFunEvals: %i  Bestf: %g   Bestfit: %g   Pen: %g   CPUTime: %f    Var: %g \n',nfuneval,fbest,fbest - reg_alpha*objReg,objReg, cputime-cpu_time,var(Refset_values_penalty));
    % else
%% Nov2016: correction when there is no feasible point    
     elseif isinf(fbest)
        fprintf('Initial Pop: NFunEvals: %i  NO FEASIBLE POINT FOUND  Bestf (penalized): %g      CPUTime: %f    Var: %g \n',nfuneval,Refset_values_penalty(1),cputime-cpu_time,var(Refset_values_penalty));
     else
        fprintf('Initial Pop: NFunEvals: %i  Bestf: %g      CPUTime: %f    Var: %g \n',nfuneval,fbest,cputime-cpu_time,var(Refset_values_penalty));
     end
       
    %end
end %% iterprint==1

Results.f=[Results.f fbest];
Results.x=[Results.x;xbest];
Results.time=[Results.time cputime-cpu_time];
Results.neval=[Results.neval nfuneval];

%Sort the Refset
[Refset_values_penalty,I]=sort(Refset_values_penalty);
Refset_values=Refset_values(I);
penalty=penalty(I);
Refset=Refset(I,:);
Refset_nlc=Refset_nlc(I,:);

if combin==1, nrand=nvar; elseif combin==2, nrand=1; end

final_ref=0;
use_bestx=0;

index1=ncomb(:,1);
index2=ncomb(:,2);
index=[index1;index2];


diff_index=index2-index1;

lb_p=1;         %Minimum ppp
st_p=0.75;      %Defines maximum ppp. max(ppp)= lb_p+ st_p
%Example: lb_p=1.25 and st_p=0.5 varies ppp from
%1.25 to 1.75

ppp=st_p*(diff_index-1)/(dim_refset-2)+lb_p;
hyper_x_L=repmat(x_L,MaxSubSet,1);
hyper_x_U=repmat(x_U,MaxSubSet,1);

refset_change=zeros(1,dim_refset);

while (not(fin))
    child=zeros(MaxSubSet2,nvar);
    child_values=zeros(MaxSubSet2,1);
    child_values_penalty=zeros(MaxSubSet2,1);
    child_penalty=zeros(MaxSubSet2,1);
    child_nlc=zeros(MaxSubSet2,nconst);
    
    child_parent_index=zeros(MaxSubSet2,1);
    counter=1;
    continuar=0;
    while ~continuar
        %Sort Refset
        %This is mandatory because the variable index is also sorted
        [Refset_values_penalty,I]=sort(Refset_values_penalty);
        Refset_values=Refset_values(I);
        penalty=penalty(I);
        Refset=Refset(I,:);
        Refset_nlc=Refset_nlc(I,:);
        
        refset_change=refset_change(I);
        
        parents_index1=Refset(index1,:);
        parents_index2=Refset(index2,:);
        parents_index1_values_penalty=Refset_values_penalty(index1);
        parents_index2_values_penalty=Refset_values_penalty(index2);
        
        denom22=max(abs(parents_index1),abs(parents_index2));
        
        denom22(find(~denom22))=1;
        
        AAA=abs((parents_index1-parents_index2)./denom22);
        BBB=find(all(AAA'<1e-3));
        if ~isempty(BBB)
            %Check the index of worst element among the twins of index2
            index_refset_out=max(index2(BBB));
            include=0;
            while ~include
                %Generate a new random solution
                new_refset_member=rand(1,nvar).*(xu_log-xl_log)+xl_log;
				
				if log_level < 2
                new_refset_member(log_var)=exp(new_refset_member(log_var));
                end
				
                [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(new_refset_member,x_L,x_U,fobj,nconst,c_L,c_U,...
                    tolc,weight,int_var,bin_var,log_level,log_var,varargin{:});
                nfuneval=nfuneval+1;
                if include
                    Refset(index_refset_out,:)=x;
                    Refset_values(index_refset_out)=val;
                    Refset_values_penalty(index_refset_out)=val_penalty;
                    Refset_nlc(index_refset_out,:)=nlc;
                    penalty(index_refset_out)=pena;
                    members_update(index_refset_out)=0;
                end
            end
        else
            continuar=1;
        end
    end
    
    candidate=Refset;
    candidate_values=Refset_values;
    candidate_values_penalty=Refset_values_penalty;
    candidate_nlc=Refset_nlc;
    candidate_penalty=penalty;
    candidate_update=zeros(1,dim_refset);
    
    members_update=ones(1,dim_refset);
    
    factor=repmat(ppp,1,nvar).*(parents_index2-parents_index1)/1.5;
    
    v1=parents_index1-factor;
    v2=parents_index2-factor;
    v3=2*parents_index2-parents_index1-factor;
    
    aaa=find(v1<hyper_x_L);
    if ~isempty(aaa)
        rand_aaa=rand(1,length(aaa));
        AAA=find(rand_aaa>prob_bound);
        aaa=aaa(AAA);
        v1(aaa)=hyper_x_L(aaa);
    end
    
    bbb=find(v1>hyper_x_U);
    if ~isempty(bbb)
        rand_bbb=rand(1,length(bbb));
        BBB=find(rand_bbb>prob_bound);
        bbb=bbb(BBB);
        v1(bbb)=hyper_x_U(bbb);
    end
    
    ccc=find(v3<hyper_x_L);
    if ~isempty(ccc)
        rand_ccc=rand(1,length(ccc));
        CCC=find(rand_ccc>prob_bound);
        ccc=ccc(CCC);
        v3(ccc)=hyper_x_L(ccc);
    end
    
    ddd=find(v3>hyper_x_U);
    if ~isempty(ddd)
        rand_ddd=rand(1,length(ddd));
        DDD=find(rand_ddd>prob_bound);
        ddd=ddd(DDD);
        v3(ddd)=hyper_x_U(ddd);
    end
    
    if nrand>1
        new_comb1=v1+(v2-v1).*rand(MaxSubSet,nrand);
        new_comb2=v2+(v3-v2).*rand(MaxSubSet,nrand);
    else
        new_comb1=v1+(v2-v1)*rand;
        new_comb2=v2+(v3-v2)*rand;
    end
    
    new_comb=[new_comb1;new_comb2];
    
    for i=1:MaxSubSet2
        %Evaluate
        [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(new_comb(i,:),x_L,x_U,fobj,nconst,c_L,c_U,...
            tolc,weight,int_var,bin_var, log_level, log_var,varargin{:});
        nfuneval=nfuneval+1;
        
        if include
            child(counter,:)=x;
            child_values(counter)=val;
            child_values_penalty(counter)=val_penalty;
            child_penalty(counter)=pena;
            child_nlc(counter,:)=nlc;
            child_parent_index(counter)=index(i);
            
            counter=counter+1;
            
            if val_penalty<candidate_values_penalty(index(i))
                %Update candidate
                %This can be accelerated saving the index and replacing at
                %the end with the best one
                candidate(index(i),:)=x;
                candidate_values_penalty(index(i))=val_penalty;
                candidate_values(index(i))=val;
                candidate_penalty(index(i))=pena;
                candidate_nlc(index(i),:)=nlc;
                candidate_update(index(i))=1;
                members_update(index(i))=0;
            end
        end
    end
    
    
    %Check the stop criterion
    if nfuneval>= maxeval
        fin=1;
    elseif cputime-cpu_time>=maxtime
        fin=2;
    elseif stopOptimization
        fin=4;
    elseif not(isempty(vtr))
        if fbest<=vtr
            fin=3;
        end
    end
    
    child(counter:end,:)=[];
    child_values(counter:end)=[];
    child_values_penalty(counter:end)=[];
    child_penalty(counter:end)=[];
    child_parent_index(counter:end)=[];
    
    
    members_to_update=find(candidate_update==1);
    for i=1:length(members_to_update);
        [vector,vector_value,vector_penalty,vector_value_penalty,vector_nlc,new_child,new_child_value,...
            new_child_penalty,new_child_value_penalty,new_child_nlc,nfuneval]=ssm_beyond(Refset(members_to_update(i),:),...
            candidate(members_to_update(i),:),candidate_values_penalty(members_to_update(i)),fobj,nrand,tolc,weight,...
            x_L,x_U,c_L,c_U,nconst,int_var,bin_var,nfuneval,prob_bound,log_level,log_var,varargin{:});
        
        if not(isempty(vector))
            Refset(members_to_update(i),:)=vector;
            Refset_values(members_to_update(i))=vector_value;
            Refset_values_penalty(members_to_update(i))=vector_value_penalty;
            penalty(members_to_update(i))=vector_penalty;
            Refset_nlc(members_to_update(i),:)=vector_nlc;
        else
            Refset(members_to_update(i),:)=candidate(members_to_update(i),:);
            Refset_values(members_to_update(i))=candidate_values(members_to_update(i));
            Refset_values_penalty(members_to_update(i))=candidate_values_penalty(members_to_update(i));
            penalty(members_to_update(i))=candidate_penalty(members_to_update(i));
            Refset_nlc(members_to_update(i),:)=candidate_nlc(members_to_update(i),:);
        end
    end
    
    fbest_lastiter=fbest;
    
    ggg=find(penalty<=tolc);
    if not(isempty(ggg))
        [fbest_new,hhh]=min(Refset_values_penalty(ggg));
        if fbest_new<fbest
            xbest=Refset(ggg(hhh),:);
            fbest=fbest_new;
            use_bestx=1;
            if inter_save
                Results.fbest=fbest;
                Results.xbest=xbest;
                Results.numeval=nfuneval;
                Results.end_crit=fin;
                Results.cpu_time=cpu_time;
                Results.Refset.x=Refset;
                Results.Refset.f=Refset_values;
                Results.Refset.fpen=Refset_values_penalty;
                Results.Refset.const=Refset_nlc;
                Results.Refset.penalty=penalty;
                save ess_report problem opts Results
            end
        end
    end
    
    iter=iter+1;
    n_minimo=n_minimo+1;
    
    if iterprint==1
        if reg_flag
            % calculate the regularization term to substract from the objective.
            amigo_privstruct.theta = xbest;
            amigo_privstruct.theta = AMIGO_scale_theta_back(amigo_privstruct.theta, amigo_inputs);
            [objReg] = AMIGO_PEcostreg(amigo_privstruct.theta,amigo_inputs,amigo_privstruct);
            reg_alpha = amigo_inputs.nlpsol.regularization.alpha;
            fprintf('Iteration: %i NFunEvals: %i  Bestf: %g   Bestfit: %g   Pen: %g   CPUTime: %f    Var: %g \n',iter,nfuneval,fbest,fbest - reg_alpha*objReg,objReg, cputime-cpu_time,var(Refset_values_penalty));
           %else
           %% correction Nov-2016 - check if there are no feasible points
            elseif isinf(fbest)
            fprintf('Iteration: %i NFunEvals: %i  NO FEASIBLE POINT FOUND  Best Refset value (penalized): %g      CPUTime: %f    Var: %g \n',iter,nfuneval,Refset_values_penalty(1),cputime-cpu_time,var(Refset_values_penalty));
            else
             fprintf('Iteration: %i NFunEvals: %i  Bestf: %g      CPUTime: %f    Var: %g \n',iter,nfuneval,fbest,cputime-cpu_time,var(Refset_values_penalty));
            end
           %% fprintf('Iteration: %i NFunEvals: %i  Bestf: %g      CPUTime: %f    Var: %g \n',iter,nfuneval,fbest,cputime-cpu_time,var(Refset_values_penalty));
        %end
    end
    
    ddd=find(members_update);
    eee=find(~members_update);
    refset_change(ddd)=refset_change(ddd)+1;
    refset_change(eee)=0;
    fff=find(refset_change>=20);
    if ~isempty(fff)
        to_replace=numel(fff);
        replaced=0;
        while replaced<to_replace
		% BUG FIXED JUNE-2014 - JRB
		%            newsol=rand(1,nvar).*(xu_log-xl_log)+xl_log;
            newsol=rand(1,nvar).*(xu_log-xl_log)+xl_log;
			if log_level < 2
			newsol(log_var)=exp(newsol(log_var)); 
			end
            
            [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(newsol,x_L,x_U,fobj,nconst,c_L,c_U,...
                tolc,weight,int_var,bin_var,log_level,log_var,varargin{:});
            nfuneval=nfuneval+1;
            if include
                Refset(fff(replaced+1),:)=x;
                Refset_values(fff(replaced+1))=val;
                Refset_values_penalty(fff(replaced+1))=val_penalty;
                penalty(fff(replaced+1))=pena;
                Refset_nlc(fff(replaced+1),:)=nlc;
                refset_change(fff(replaced+1))=0;
                replaced=replaced+1;
            end
        end
    end
    
    if fbest<fbest_lastiter %((fbest<fbest_lastiter)/fbest_lastiter)>tolc
        Results.f=[Results.f fbest];
        Results.x=[Results.x; xbest];
        Results.time=[Results.time cputime-cpu_time];
        Results.neval=[Results.neval nfuneval];
        nreset=0;
        if plot_results==1
            stairs(Results.time,Results.f)
            xlabel('CPU Time (s)');
            ylabel('Objective Function Value (-)');
            title(strcat('Best f:  ',num2str(fbest)))
            drawnow
        end
    else
        if n_stuck
            nreset=nreset+1;
            if nreset==n_stuck
                fin=5;
            end
        end
    end
    
    if local_solver & local_bestx & n_minimo>=n_critico
        if use_bestx
            ess_local_filters
            n_minimo=0;
        end
    else
        if (local_solver & n_minimo>=n_critico)
            ess_local_filters
            n_minimo=0;
        end
    end
    fbest_lastiter=fbest;
    
    %Check the stopping criterion
    if nfuneval>= maxeval
        fin=1;
    elseif cputime-cpu_time>=maxtime
        fin=2;
    elseif stopOptimization
        fin=4;
    elseif not(isempty(vtr))
        if fbest<vtr
            fin=3;
        end
    end
    
    if fin
        fprintf('\n');
        if local_finish & fin<3
            final_ref=1;
            if isempty(strmatch(local_solver,local_finish))
                %ssm_aux_local(problem.f,x_L,x_U,c_L,c_U,neq,local_finish,nvar,varargin{:});
                if isFjacDefined
                    ssm_aux_local(problem.f,problem.fjac,x_L,x_U,c_L,c_U,neq,local_finish,nvar,varargin{:});
                else
                    ssm_aux_local(problem.f,[],x_L,x_U,c_L,c_U,neq,local_finish,nvar,varargin{:});
                end
            end
            local_tol=3;
            
            if not(isempty(xbest))
                x0=xbest;
                f0=fbest;
            else
                x0=Refset(1,:);
                f0=Refset_values_penalty(1);
            end
            if iterprint
                fprintf('Final local refinement with: %s \n', upper(local_finish))
                fprintf('Initial point function value: %f \n',f0);
                tic
            end
            if isFjacDefined
                [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,fjac,...
                    local_finish,local_iterprint,local_tol,weight,nconst,tolc,opts.local,log_level,log_var,varargin{:});
            else
                [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,[],...
                    local_finish,local_iterprint,local_tol,weight,nconst,tolc,opts.local,log_level,log_var,varargin{:});
            end
            %      [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,...
            %    local_finish,local_iterprint,local_tol,weight,nconst,tolc,varargin{:});
            
            if iterprint
                fprintf('Local solution function value: %f \n',fval);
                fprintf('Number of function evaluations in the local search: %i \n',numeval);
                fprintf('CPU Time of the local search: %f  seconds \n\n',toc);
            end
            nfuneval=nfuneval+numeval;
            %Evaluate the local solution
            [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(x,x_L,x_U,fobj,nconst,c_L,c_U,...
                tolc,weight,int_var,bin_var,log_level,log_var,varargin{:});
            nfuneval=nfuneval+1;
            if include & pena<=tolc
                if val_penalty<fbest
                    fbest=val_penalty;
                    xbest=x;
                    if fbest<vtr, fin=3; end
                end
            end
            ssm_delete_files(local_finish,c_U);
        end
        cpu_time=cputime-cpu_time;
        Results.f=[Results.f fbest];
        Results.x=[Results.x; xbest];
        Results.neval=[Results.neval nfuneval];
        Results.time=[Results.time cpu_time];
        if iterprint %%
            
        %%% Jake Fixes double display of results and displays all results in linear scale % Start
        if isinf(fbest)
            
        fprintf('NO FEASIBLE SOLUTION FOUND   Best Refset value (penalized)\t\t%g\n', Refset_values_penalty(1));
        fprintf('Decision vector\n');
                
            if log_level > 1
                
            Refset_dummy = Refset(1,:)';
            
            Refset_dummy(log_var) = exp(Refset_dummy(log_var));
                
            fprintf('\t%g\n', Refset_dummy);
             
            else
                
            fprintf('\t%g\n', Refset(1,:)');
        
            end
                      
        else 
            
            switch fin
                case 1
                    disp('Maximum number of function evaluations achieved')
                case 2
                    disp('Maximum computation time achieved')
                case 3
                    disp('Desired function value achieved')
                case 4
                    disp('Optimization stopped by the user')
                case 5
                    disp('Maximum number of iterations without significant improvement')
            end %%switch
            fprintf('Best solution value\t\t%g\n', fbest);
            fprintf('Decision vector\n');
            
            if log_level > 1
                
            xbest_dummy = xbest';
            
            xbest_dummy(log_var) = exp(xbest_dummy(log_var));
                
            fprintf('\t%g\n', xbest_dummy);
             
            else
                
            fprintf('\t%g\n', xbest');        
            
            end
            
            
        end        
        %%% Jake Fixes double display of results and displays all results in linear scale % End
        
%% Correction Nov2016 % Jake 2018 see above
           %if iterprint==1
%                 if isinf(fbest)
%                     fprintf('NO FEASIBLE SOLUTION FOUND   Best Refset value (penalized)\t\t%g\n', Refset_values_penalty(1));
%                     fprintf('Decision vector\n');
%                     fprintf('\t%g\n', Refset(1,:)');
%                 else
%                     fprintf('Best solution value\t\t%g\n', fbest);
%                     fprintf('Decision vector\n');
%                     fprintf('\t%g\n', xbest');
%                 end
            %end
           
            fprintf('CPU time\t\t%g\n', cpu_time);
            fprintf('Number of function evaluations\t\t%g\n', nfuneval);
        end %% if iterprint
        fclose all;
        Results.fbest=fbest;
        Results.xbest=xbest;
        Results.numeval=nfuneval;
        Results.end_crit=fin;
        Results.cpu_time=cpu_time;
        Results.Refset.x=Refset;
        Results.Refset.f=Refset_values;
        Results.Refset.fpen=Refset_values_penalty;
        Results.Refset.const=Refset_nlc;
        Results.Refset.penalty=penalty;
        if plot_results==1 | plot_results==2
            stairs(Results.time,Results.f)
            xlabel('CPU Time (s)');
            ylabel('Objective Function Value (-)');
            title(strcat('Best f:  ',num2str(fbest)))
        end
        try
            save ess_report problem opts Results
        end
        %Delete aux files
        if local_solver & isempty(strmatch(local_solver,local_finish));
%             ssm_delete_files(local_solver,c_U);
        end
        %     keyboard
        %%% Jake Output results in linear scale. % Start
        if log_level > 1
            
        pos = zeros([1, length(log_var)]);   
        
        pos(log_var) = 1;
        
        pos = repmat(pos, size(Results.x, 1), 1) == 1;
        
        Results.x(pos) = exp(Results.x(pos));
            
        Results.xbest(log_var) = exp(Results.xbest(log_var));
        
        pos = zeros([1, length(log_var)]);   
        
        pos(log_var) = 1;
        
        pos = repmat(pos, size(Results.Refset.x, 1), 1) == 1;
        
        Results.Refset.x(pos) = exp(Results.Refset.x(pos));
        %%% Jake Output results in linear scale. % End

        end
        
        return
        
    end
end

%\-----------------------------------------------------/
% Definition of objective for fminsearchbnd
function fp = fminsobj(x,fun,c_L,c_U,weight,tolc)

[f,c] = feval(fun,x);
pen=ssm_penalty_function(x,c,c_L,c_U,tolc);

fp=f+weight*pen;

%n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/

