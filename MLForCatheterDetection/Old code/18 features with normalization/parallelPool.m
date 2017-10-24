tic;
delete(gcp('nocreate'))
poolSize = 3;
profileName = 'local';
pool = parpool(profileName,poolSize);
poolInitializationTime = toc;
job = batch('netTrain', 'Pool', poolSize);
sprintf('Pool initialization time is %.2f seconds', poolInitializationTime)


