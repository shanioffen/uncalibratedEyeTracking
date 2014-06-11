% so - July 12 2011, modified from Justin's testExperiment

function [myscreen task stimulus] = EyetrackSimple(doParam,movieChoice)
global stimulus;

if(nargin<2)
    movieChoice = 2;
end

stimulus.movieChoice = movieChoice;

whichEye = input('which eye? [R or L]','s');
stimulus.whichEye = whichEye;


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

% set important variables
setMovieDetails;

mglSetParam('movieMode',1,1);
segLength = stimulus.movieLength;

task{1}.collectEyeData = 1;
task{1}.waitForBacktick = 0;
task{1}.seglen = [segLength];
task{1}.numBlocks = 1;
task{1}.numTrials = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'RUN THE TASK PROGRAM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the main task

[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@updateScreenCallback);

  
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

saveSubjData(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

myscreen = mytrack(myscreen);

if(stimulus.startR)
 mglEyelinkRecordingStart(myscreen.eyetracker.data);
 stimulus.m = mglMovie(stimulus.movieFile);
 mglMovie(stimulus.m,'play');   
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      individualized  SUBFUNCTIONS                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setMovieDetails;
global stimulus;

movieDir = {'/Users/shani/Movies/'};

movieList = {'CDwaka_30pct_20crop_10s.mov','CWandCCWvogueWcount.mov','fantasiaCCW_wCD.mov'};
movieChoice = stimulus.movieChoice;
movie = movieList{movieChoice};

movieLengthList ={400,400,400};
movieLength = movieLengthList{movieChoice};

movieFile = [movieDir{1} movie];
stimulus.movieFile = movieFile;
stimulus.movie = movie;
stimulus.movieLength = movieLength;   
   
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function saveSubjData(myscreen)
  global stimulus
  
    prompt = {...
            '1: movieType',...
            '2: subjectCode',...
            '3: EDF filename',...
            '4: Todays date and time'...
            '5: Control or Test?'...
           };
  
  defaultMovie =  ...
  [stripext(stimulus.movie) '_' num2str(stimulus.movieLength)];
  
  edfName = [num2str(myscreen.eyetracker.datafilename) '.edf'];
                   
  tempname = 'Set subject information to save with data';
  numlines = 2; 
  defaultanswer = {...
                   defaultMovie,...
                   num2str(round(1000000*rand)),...
                   edfName,...
                   num2str(clock),...
                   'Test',...
                  };
  
  options.Resize = 'on';
  options.WindowStyle = 'normal';
  options.Interpreter = 'none';
  
 
  
  SubjInfo.answer = defaultanswer;
  disp(['Subj code num is ' num2str(defaultanswer{2})])
  

  
  
  SubjInfo.movieType = SubjInfo.answer{1};
  SubjInfo.subjectCode = SubjInfo.answer{2};
  SubjInfo.STIMfilename = [SubjInfo.answer{2} '.mat'];
  SubjInfo.visitNum = SubjInfo.answer{3};
  SubjInfo.todaysDate = SubjInfo.answer{4};
  SubjInfo.CorT = [SubjInfo.answer{5}];
  stimulus.SubjInfo = SubjInfo;
  



  
  stimulus.dataDir = '/Users/shani/data/stimulusVariable/';
  cd(stimulus.dataDir);
  save(SubjInfo.STIMfilename, 'stimulus');
  disp('done saving stim file');

