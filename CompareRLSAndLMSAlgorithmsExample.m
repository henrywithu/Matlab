%% Compare RLS and LMS Algorithms
% Equalize a QAM signal passed through a frequency-selective fading
% channel using RLS and LMS algorithms. Compare the performance of the
% two algorithms.
%%
% Specify the modulation order. Generate the corresponding QAM reference
% constellation.
M = 16; 
sigConst = qammod(0:M-1,M,'UnitAveragePower',true);
%%
% Create a frequency-selective static channel having three taps.
rchan = comm.RayleighChannel('SampleRate',1000, ...
    'PathDelays',[0 1e-3 2e-3],'AveragePathGains',[0 -3 -6], ...
    'MaximumDopplerShift',0, ...
    'RandomStream','mt19937ar with seed','Seed',73);

%% RLS Equalizer
% Create an RLS equalizer object.
eqrls = lineareq(6,rls(0.99,0.1)); 
eqrls.SigConst = sigConst; 
eqrls.ResetBeforeFiltering = 0; 
%%
% Generate and QAM modulate a random training sequence. Pass the
% sequence through the Rayleigh fading channel. Pass the received signal
% and the training signal through the equalizer to set the equalizer tap
% weights.
trainData = randi([0 M-1],200,1);
trainSig = qammod(trainData,M,'UnitAveragePower',true);
rxSig = rchan(trainSig);
[~,~,errorSig] = equalize(eqrls,rxSig,trainSig);
%%
% Plot the magnitude of the error estimate. 
plot(abs(errorSig))
title('Error Estimate, RLS Equalizer')
xlabel('Symbols')
ylabel('Amplitude')
%%
% The error is nearly eliminated within 200 symbols.
%%
% Transmit a QAM signal through a frequency-selective channel. Equalize the
% received signal using the previously 'trained' RLS equalizer. Measure the
% time required to execute the processing loop.
tic
for k = 1:20
   data = randi([0 M-1],1000,1); % Random message
   txSig = qammod(data,M,'UnitAveragePower',true);

   % Introduce channel distortion.
   rxSig = rchan(txSig);

   % Equalize the received signal.
   eqSig = equalize(eqrls,rxSig);

end
rlstime = toc;
%%
% Plot the constellation diagram of the received and equalized signals. 
h = scatterplot(rxSig,1,0,'c.');
hold on
scatterplot(eqSig,1,0,'b.',h)
legend('Received Signal','Equalized Signal')
title('RLS Equalizer')
hold off
%%
% The equalizer removed the effects of the fading channel.
%% LMS Equalizer
% Repeat the equalization process with an LMS equalizer. Create an LMS
% equalizer object.
eqlms = lineareq(6,lms(0.03)); 
eqlms.SigConst = sigConst; 
eqlms.ResetBeforeFiltering = 0; 
%%
% Train the LMS equalizer.
trainData = randi([0 M-1],1000,1);
trainSig = qammod(trainData,M,'UnitAveragePower',true);
rxSig = rchan(trainSig);
[~,~,errorSig] = equalize(eqlms,rxSig,trainSig);
%%
% Plot the magnitude of the error estimate. 
plot(abs(errorSig))
title('Error Estimate, LMS Equalizer')
xlabel('Symbols')
ylabel('Amplitude')
%%
% Training the LMS equalizer requires 1000 symbols.
%%
% Transmit a QAM signal through the same frequency-selective channel.
% Equalize the received signal using the previously 'trained' LMS
% equalizer. Measure the time required to execute the processing loop.
tic
for k = 1:20
   data = randi([0 M-1],1000,1); % Random message
   txSig = qammod(data,M,'UnitAveragePower',true);

   % Introduce channel distortion.
   rxSig = rchan(txSig);

   % Equalize the received signal.
   eqSig = equalize(eqlms,rxSig);

end
lmstime = toc;
%%
% Plot the constellation diagram of the received and equalized signals. 
h = scatterplot(rxSig,1,0,'c.');
hold on
scatterplot(eqSig,1,0,'b.',h)
legend('Received Signal','Equalized Signal')
title('LMS Equalizer')
%%
% The equalizer removes the effects of the fading channel.
%%
% Compare the loop execution time for the two equalizer algorithms.
[rlstime lmstime]
%%
% The LMS algorithm is more computationally efficient as it took 50% of the
% time to execute the processing loop. However, the training sequence
% required by the LMS algorithm is 5 times longer.