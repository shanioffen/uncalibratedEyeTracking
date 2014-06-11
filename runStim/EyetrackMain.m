function EyetrackMain
  global stimulus;
  


% set important variables
setMovieDetails;


% call the task script that shows the movie
  %*%*%*%*%*%%*%*%*%*%*%*%*%%*%*%*%%%*%%*%*%*%*%%%%*&&%&*^*%&% 
[myscreen, task, stimulus] = EyetrackSubsidiary(0);
   %*%*%*%*%*%%*%*%*%*%*%*%*%%*%*%*%%%*%%*%*%*%*%%%%*&&%&*^*%&% 
   
% make any necessary adjustments, and save  
saveSubjData(myscreen);
   
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      individualized  SUBFUNCTIONS                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setMovieDetails;
global stimulus;

movieDir = {'/Users/shani/Movies/'};

movieList = {'shakira.mov','CDwaka_30pct_20crop_10s.mov'};
movieChoice = 2;
movie = movieList{movieChoice};

movieLengthList ={360,400};
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
  
  % Answer = ...
  %      inputdlg(prompt,tempname,numlines,defaultanswer,options);

  
  %  if isempty(answer);
  %   display('Exited without completing the form; using defaults');
  %   SubjInfo.answer = defaultanswer;
  % else
  %   SubjInfo.answer = answer;    
  % end
  
  SubjInfo.answer = defaultanswer;
  disp(['Subj code num is ' num2str(defaultanswer{2})])
  
  whichEye = input('which eye? [R or L]','s');
  subjInfo.whichEye = whichEye;

  
  
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
