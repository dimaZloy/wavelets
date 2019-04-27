#WAVELET  1D Wavelet transform with optional singificance testing
#
#   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
#
# INPUTS:
#
#    Y = the time series of length N.
#    DT = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)^2.
#    Its units are sigma^2 (the time series variance).
#
#
# OPTIONAL INPUTS:
# 
# *** Note *** setting any of the following to -1 will cause the default
#               value to be used.
#
#    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    DJ = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
#        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
#
#    MOTHER = the mother wavelet function.
#             The choices are 'MORLET', 'PAUL', or 'DOG'
#
#    PARAM = the mother wavelet parameter.
#            For 'MORLET' this is k0 (wavenumber), default is 6.
#            For 'PAUL' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#
# OPTIONAL OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
#            where J1+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot COI lines on a contour plot by doing:
#
#              contour(time,log(period),log(power))
#              plot(time,log(coi),'k')
#
#----------------------------------------------------------------------------
#   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made. This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------

#include("wave_bases.jl");

function wavelet(Y,dt,pad=0,dj=-1,s0=-1,J1=-1,mother=-1,param=-1)

	wave = period = scale = coi = fourier_factor = 0;

	n1 = length(Y);

	if (s0 == -1)
		s0=2.0*dt;
	end
	if (dj == -1)
		dj = 1.0/4.0;
	end
	if (J1 == -1)
		# J1= fix((log(n1*dt/s0)/log(2))/dj);
		J1= round( Int,  (log(n1*dt/s0)/log(2))/dj ) ;
	end
	if (mother == -1)
		mother = "MORLET";
	end

	#....construct time series to analyze, pad if necessary
	# x(1:n1) = Y - mean(Y);
	x = zeros(Float64,n1);
	x = Y .-mean(Y);
	
	
	if (pad == 1)
		
		# base2 = fix(log(n1)/log(2) + 0.4999);   # power of 2 nearest to N
		base2 = round( Int, log(n1)/log(2) + 0.4999 );
		
		#x = [x,zeros(1,2^(base2+1)-n1)];
		x = vcat(x, zeros(2^(base2+1)-n1,1));
		
	end
	
	n = length(x);

	#....construct wavenumber array used in transform [Eqn(5)]
	#k = [1:fix(n/2)];
	#k = k.*((2.*pi)/(n*dt));
	#k = [0., k, -k(fix((n-1)/2):-1:1)];
	
	###k = LinRange(1, round(Int, n/2), round(Int, n/2) ];
	k = LinRange(1: n/2);
	k = k.*( (2.0*pi)/(n*dt) );
	kInv  = -reverse(k);
	kInv =  kInv[2: round(Int, (n-1)/2)];
	k = vcat(0.0,k, kInv);
	

	#....compute FFT of the (padded) time series
	f = fft(x);    # [Eqn(3)]

	#....construct SCALE array & empty PERIOD & WAVE arrays
	J1i = round(Int,J1);	
	a = LinRange(0,J1i,J1i+1).*dj;	
	scale = ones(length(a),1);
	for i=1:length(a)
		scale[i] = s0*2.0^a[i];		
	end
	
	period = scale;
		
		
		
	wave = zeros( round(Int, J1+1),n );  # define the wavelet array
	#wave = zeros(Float64, J1+1);	
	wave = wave .+ wave.*im;  # make it complex

	# loop through all scales and compute transform
	for a1 = 1: round(Int, J1+1)
		daughter,fourier_factor,coi,dofmin =wave_bases(mother,k,scale[a1],param);	
		wave[a1,:] = ifft(f.*daughter);  # wavelet transform[Eqn(4)]
	end

	period = fourier_factor.*scale;
	coi = coi*dt*[1E-5,1:((n1+1)/2-1), reverse((1:(n1/2-1))), 1E-5];  # COI [Sec.3g]
	wave = wave[:,1:n1];  # get rid of padding before returning

	return wave,period,scale,coi; 

end # of code

#WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
#
#  [DAUGHTER,FOURIER_FACTOR,COI,DOFMIN] = ...
#      wave_bases(MOTHER,K,SCALE,PARAM);
#
#   Computes the wavelet function as a function of Fourier frequency,
#   used for the wavelet transform in Fourier space.
#   (This program is called automatically by WAVELET)
#
# INPUTS:
#
#    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
#    K = a vector, the Fourier frequencies at which to calculate the wavelet
#    SCALE = a number, the wavelet scale
#    PARAM = the nondimensional parameter for the wavelet function
#
# OUTPUTS:
#
#    DAUGHTER = a vector, the wavelet function
#    FOURIER_FACTOR = the ratio of Fourier period to scale
#    COI = a number, the cone-of-influence size at the scale
#    DOFMIN = a number, degrees of freedom for each point in the wavelet power
#             (either 2 for Morlet and Paul, or 1 for the DOG)
#
#----------------------------------------------------------------------------
#   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
#   University of Colorado, Program in Atmospheric and Oceanic Sciences.
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made.  This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#----------------------------------------------------------------------------
function wave_bases(mother,k,scale,param);

	mother = uppercase(mother);
	n = length(k);

	if (mother == "MORLET")  #-----------------------------------  Morlet
		
		if (param == -1)
			param = 6.;
		end
		
		k0 = param;
		kpos = zeros(length(k),1);
		id = findall(x->x>0.0,k);
		for i=1:length(id)
			kpos[id[i]] = 1;	
		end
		
		expnt = -(scale.*k .- k0).^2/2;
		expnt = expnt.*kpos;
			
		norm = sqrt(scale*k[2])*(pi^(-0.25))*sqrt(n);    # total energy=N   [Eqn(7)]
		daughter = norm.*exp.(expnt);
		daughter = daughter.*kpos;     # Heaviside step function
		fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); # Scale-->Fourier [Sec.3h]
		coi = fourier_factor/sqrt(2);                  # Cone-of-influence [Sec.3g]
		dofmin = 2;                                    # Degrees of freedom
	elseif (mother == "PAUL")  #--------------------------------  Paul
		#if (param == -1), param = 4.;, end
		#m = param;
		#expnt = -(scale.*k).*(k > 0.);
		#norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
		#daughter = norm*((scale.*k).^m).*exp(expnt);
		#daughter = daughter.*(k > 0.);     # Heaviside step function
		#fourier_factor = 4*pi/(2*m+1);
		#coi = fourier_factor*sqrt(2);
		#dofmin = 2;
		error("PAUL is not implemented")
	elseif (mother == "DOG")  #--------------------------------  DOG
		#if (param == -1), param = 2.;, end
		#m = param;
		#expnt = -(scale.*k).^2 ./ 2.0;
		#norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
		#daughter = -norm*(i^m)*((scale.*k).^m).*exp(expnt);
		#fourier_factor = 2*pi*sqrt(2./(2*m+1));
		#coi = fourier_factor/sqrt(2);
		#dofmin = 1;
		error("DOG is not implemented")
	else
		error("Mother must be one of MORLET,PAUL,DOG")
	end

	return daughter,fourier_factor,coi,dofmin; 

end # of code

function computeWaveletContourMatrix(time::Array{Float64,1}, signal::Array{Float64,1},dt::Float64, UInf::Float64,D::Float64 )
	
	pad = 1;
	dj  = 0.2;
	so = 2.0*dt;
	j1 = -1;
	lag1 = 0.72;
	mother = "Morlet";
	
	wave, period, scale, coi = wavelet(signal, dt, pad, dj,so, j1,mother);
	
	power = (abs.(wave)).^2;

	nx = length(signal);
	ny = length(period);

	tc = zeros(Float64,ny,nx);
	yc = zeros(Float64,ny,nx);

	for j = 1:ny
		tc[j,1:end] = time*UInf/D;
	end

	for j =1:nx
		for k = 1:ny
			yc[k,j] = 1.0/period[k];
		end
	end

	return tc, yc, power ; 

end
