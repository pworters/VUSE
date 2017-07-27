%
%   function [auc,vuseamp,flipScheme,flipOpt] = 
%                               make_vuse(np, TR, T1, T2, nramp, plotfigs)
%
%	Function modified from Brian's SSFP code.  This code generates a
%	flipScheme that maintains a constant echo amplitude over the SSFP
%	train.  The echo amplitude target is recursively determined and
%   maximised for a certain T1, T2, TR and pulse block length.
%
%	INPUT:
%       np         = number of acq pulses (or pulse block length)
% 		T1,T2,TR all in milliseconds.
%       nramp      = number of Kaiser pulses (not included in np)
%       plotfigs   = plot RF and echoes
%
%	OUTPUTS:
%       flipScheme = array of RF pulses required for the constant echo
%       vuseamp    = targfet echo amplitude
%       flipOpt    = optimal flip angle based on T1/T2
%       auc        =  area under the echo-time curve 
%
%	NOTES:
%       Calculated for on-resonant case, and TE = TR/2 always.
%
%   EXAMPLE:
%       [auc,vamp,vflips,flipopt] = make_vuse(64,3.8,907,50,5,1);
%
%   written by Pauline Wong Worters, May 2009
%   (c) Board of Trustees, Leland Stanford Junior University
%

function [auc,vuseamp,flipScheme,flipOpt] = make_vuse(np, TR, T1, T2, ...
                                                      nramp, plotfigs)

    M0 = [0 0 1]';	    % Default initial magnetization for transient resp.
    TE = TR/2;          % always for these bSSFP simulations
    rfdphase = 180;     % RF phase cycling
    kaiserN = nramp;    % Use a Kaiser ramp starting condition

    %% Calculate optimal flip angle
    flipOpt = acos( (T1/T2-1)/(T1/T2+1) )*180/pi;

    %% Calculate steady state amplitude and use this as starting point
    etr1 = exp(-TR/T1);
    etr2 = exp(-TR/T2);
    alphar = flipOpt*pi/180; % in radians
    Mtarget = sqrt(etr2)*(1-etr1)*sin(alphar) / ...
                (1-(etr1-etr2)*cos(alphar)-etr1*etr2);

    %disp(sprintf('\nSteady state for const %.0f deg = %.4f \n', ...
    %    alpha, Mtarget));

    % These operators do not change so put them outside loop (not doing
    % changing TRs or any of that fancy stuff.  just variable flip.
    ete1 = exp(-TE/T1);
    ete2 = exp(-TE/T2);
    E1 = [ete2 0 0; 0 ete2 0; 0 0 ete1];
    D1 = [0; 0; 1-ete1];
    Rcyc1 = zrot(rfdphase*TE/TR);	% Rotation due to RF cycling.
    RcRE = Rcyc1*E1;
    RcREi = inv(RcRE);
    if RcREi(2,3)~=0,
        disp('HMMMMMMMMM'); % this shouldn't be!
        return
    end

    stepsize = 0.0005;
    maxAmpReached = 0;
    loopC = -1;
    Mupper = 1.0;       % upper bound to interval (which is infeasible)
    Mlower = Mtarget;   % lower bound (which is feasible)

    while loopC < 50,

        loopC = loopC+1;
        clear flipScheme Mplot

        % Init plotting array
        M = M0;
        Mplot=M;

        % Init the first angles based on target amplitude - approximately
        alpha = asin(Mtarget)*180/pi;
        alphaR = repmat(2*alpha,1,kaiserN);
        alphaR = kaiser_flips(alphaR,kaiserN)';
        % Do the ramp simulations
        for n = 1:kaiserN,
            Rx = xrot(alphaR(n));
            A = RcRE*Rx*RcRE;
            C = RcRE*Rx*D1+D1;
            M  = A*M + C;
            Mplot = [Mplot M];
        end

        flipScheme = alphaR;
        notfeasible = 0;    % reset condition

        % Let's try and calculate the flip scheme
        for n = 1:np,

            %% calculate alpha based on Mtarget
            R3 = RcRE*M + D1;
            R2M = RcREi(2,1)*Mtarget;
            F = R2M; G = R3(2); H = R3(3);
            fgh2 = (2*G*H - sqrt( (2*G*H)^2 - 4*(F^2-H^2)*(F^2-G^2) ) )...
                    / (2*(F^2-H^2));
            alpha = atan( fgh2 ) * 180 / pi;
            % atan return angle b/w -pi/2 to pi/2, so need to make this +ve
            if alpha < 0,
                alpha = 180 + alpha;
            end

            if ~isreal(alpha),
                %disp(sprintf('ERROR: Alpha not real anymore!! ( N = %d )',n));
                notfeasible = 1;
                break
            end

            %% simulation
            Rx = xrot(alpha);
            A = RcRE*Rx*RcRE;
            C = RcRE*Rx*D1+D1;
            M  = A*M + C;
            Mplot = [Mplot M];

            flipScheme = [flipScheme alpha];
        end

        % Set new target based on whether Mtarget was feasible or not
        if notfeasible == 1,
            % if not feasible, then store Mtarget as upper bound
            Mupper =  Mtarget;
        else
            % if feasible, store flips and M for plotting later
            flipSchemeTmp = flipScheme;
            MplotTmp = Mplot;

            Mlower = Mtarget;   % set feasible Mtarget as lower bound
            if abs(Mupper-Mlower) < stepsize,    % stop if close enough
                disp(sprintf('*** Finished calc at %d iterations',loopC));
                break
            end
        end
        %disp(sprintf('lower %.3f upper %.3f',Mlower,Mupper));
        Mtarget = mean([Mlower Mupper]);

    end

    flipScheme(1:kaiserN+np) = flipSchemeTmp(1:kaiserN+np); 
    Mplot = MplotTmp;
    vuseamp = Mtarget - stepsize;

    if (plotfigs == 1)
        figure;
        plot(flipScheme,'*-');
        title(sprintf('Target amp = %.3f',vuseamp));
        xlabel('RF number');
        ylabel('Flip angle / ^\circ');

        figure;
        x = 1:numel(Mplot(1,:));
        plot(x,Mplot(1,:),'-',x,Mplot(3,:),':');
        title('Magnetization vs Excitation');
        legend('M_{transverse}','M_{longtidunal}');
        xlabel('Echo number');
        ylabel('Magnetization / M_0');
    end

    %% AUC
    auc = sum(Mplot(1,kaiserN+2:kaiserN+np+1)) / np;

    %% Calculate a flip_scale factor for SAR calcs
    flip_scale = [0.0 0.0];
    max_flip = flipScheme(kaiserN+np);
    flip_scale(1) = sqrt(sum(flipScheme.*flipScheme)/...
                    (max_flip*max_flip*(numel(flipScheme))));
    flip_scale(2) = sum(flipScheme)/(max_flip*(numel(flipScheme)));



%
% Auxillary functions.
%

%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the x axis
%	by angle degrees.
function [M] = xrot(angle)
    c = cos(pi*angle/180);
    s = sin(pi*angle/180);
    M = [1 0 0; 0 c s; 0 -s c];

%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the z axis
%	by angle degrees.
function [M] = zrot(angle)
    c = cos(pi*angle/180);
    s = sin(pi*angle/180);
    M = [c s 0; -s c 0; 0 0 1];

% 	Sets up kaiser flips for SSFP catalyzation.
%	np is total number of pulses (including nramp);
function [allflips] = kaiser_flips(alpha_np,nramp)
    flips = alpha_np(:);
    np = numel(flips);
    ramp = [1:nramp]/(nramp+1);	% -- Linear
    kfilt = kaiser(nramp+1,2);
    ramp = cumsum(kfilt(1:end-1))/sum(kfilt); 	% skip last point, as
                                                % that is amplitude 1.
    flips(1:nramp) = ramp.*flips(1:nramp);
    allflips = flips;

