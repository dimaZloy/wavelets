
function simpleTestFFT()
	
	N = 2000;
	t = LinRange(0.0, 10.0*pi,N);
	ns = 1.0*rand(N);
	y = sin.(t+ns);
	
	T = 2.0*pi;
	freqPrime = 1.0/T;  # prime frequency
	
	f, p, freqPrimeFFT = computeFFT(t,y);
	f1, p1, freqPrimeFFT1 = computeFFT(t,y,1.0,1.0);
	
	println("Analytic Prime Freq [Hz]: ", freqPrime);
	println("Calculated Prime Freq [Hz]: ", freqPrimeFFT);
	println("Calculated1 Prime Freq [Hz]: ", freqPrimeFFT1);
	
	
end


function computeFFT(t,y)
	
	# slightly alternative formulation 
	
	L = length(t);
	L2 = Int( round(L/2) );
	Tspan = t[end] - t[1];
	Fs = Tspan/L2; 
	PW = fft(y);
	PW = abs.(PW/L2);
	PW = PW[1:L2].^2;
	freq = (0:L2-1)/L2/Fs;
	
	maxFreq,id = findmax(PW);
	maxFreqFFT = freq[id];
	
	return freq, PW, maxFreqFFT; 
	
end


function computeFFT(t,y, D::Float64 = 1.0, Uinf::Float64 = 1.0)
	

	N = length(t) # number of points
	N2 = Int( round(N/2) ); 
	# define time of interval, in  seconds
	T = t[end] - t[1]; #sampling period	
	p = abs.(fft(y)/ N2 );  # absolute value of the fft
	p = p[1: N2 ].^2;  ## take the power
	freq = (0: N2 - 1)/T; ## find the corresponding frequency in Hz
		
	maxP,id = findmax(p[2:end]);
	maxFreq = freq[id]; 
	
	St = maxFreq*D/Uinf;
	return freq, p, maxFreq, St;
	
end

function computeWelch(t,s)
	
	# t - time vector
	# s - signal vector
			
	N = length(t) # number of points
	
	# define time of interval, in  seconds
	T = t[end] - t[1]; #sampling period
	
	sub = 0.75;
	ov = 0.25;
	
	subset = Int(  round( sub*N) ); 
	noverlap = Int( round( ov*subset) ); 
		
	#dt = T*subset/N; 
	#Fs = T*subset/dt;
	
	#dt = T/N; 
	#Fs = T/dt;
	
	
	#pWelch = welch_pgram(s,subset, noverlap; onesided = true, nfft = nextfastfft( N), fs = Fs, window = nothing); 
	pWelch = welch_pgram(s,subset, noverlap; onesided = true, nfft = nextfastfft( N), fs = 1, window = nothing); 
	
	return pWelch;
	
end