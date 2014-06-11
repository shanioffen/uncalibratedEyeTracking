% so - July 12 2011, modified from Justin's testExperiment

function [myscreen task stimulus] = EyetrackSubsidiary(doParam)
global stimulus;

if(doParam)
  mglEditScreenParams;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen and stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen.datadir = 'Users/Shani/data';
myscreen = initScreen;
myscreen = initStimulus('stimulus',myscreen);


% set up eyetracker (may not need to do this)
myscreen.eyetracker.savedata = true;
myscreen.eyetracker.data = [1 1 1 0];
myscreen.eyetracker.collectEyeData = 1;

mglSetParam('movieMode',1,1);
segLength = stimulus.movieLength;

task{1}.collectEyeData = 1;
task{1}.waitForBacktick = 1;
task{1}.seglen = [1 segLength];
task{1}.numBlocks = 1;
task{1}.numTrials = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'RUN THE TASK PROGRAM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the main task

[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@updateScreenCallback);

tempFile = '/Users/shani/Movies/SEIN_S7_D3_E1.mov';
stimulus.ms = mglMovie(tempFile);

discard1 = input('start seinfeld? [1 or 0]');

if(discard1)
 mglMovie(stimulus.ms,'play');
end




  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
count = 0;

while (phaseNum <= length(task{1}))  && ~myscreen.userHitEsc
  % update tasks
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  myscreen = mytrack(myscreen);
elseif task.thistrial.thisseg == 2
   if(stimulus.startR)
     mglMovie(stimulus.ms,'pause');
     mglEyelinkRecordingStart(myscreen.eyetracker.data);
     stimulus.m = mglMovie(stimulus.movieFile);
     mglMovie(stimulus.m,'play');   
   end
end

discard1 = input('movie over? [1 or 0]');

if(discard1)
  mglMovie(stimulus.m,'pause');
  myscreen = endTask(myscreen,task);
  task = jumpsegment(task,inf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)
global stimulus;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       individualized  SUBFUNCTIONS                             % %



  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize eyetracker
function myscreen = mytrack(myscreen)
  global stimulus   
  

myscreen = initEyeTracker(myscreen,'Eyelink',0,1);
%mglEyelinkCMDPrintF('file_sample_data = RIGHT, GAZE, RAW, AREA, GAZERES, STATUS');
mglEyelinkOpen('100.1.1.1',0);
mglEyelinkSetup;
startR = input('start 0 or 1');
%startR=1;
stimulus.startR = startR;

