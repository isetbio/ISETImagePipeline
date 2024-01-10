function dataOut = nreConeResponse(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of cone excitations witout eye movements
%
% Syntax:
%   dataOut = nrePhotopigmentExcitationsWithNoEyeMovements(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments,
%           dataOut = nrePhotopigmentExcitationsWithNoEyeMovements()
%       it does not compute anything and simply returns a struct with the
%       defaultParams (optics and coneMosaic params) that define the neural
%       compute pipeline for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object,
%       it computes 'instancesNum' of cone photopigment excitation sequences
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%    instancesNum                   - the number of response instances to compute
%
% Optional key/value input arguments:
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are:
%                                        - 'none' (noise-free responses)
%                                        - 'random' (noisy response instances)
%                                     Default is {'random'}.
%   'rngSeed'                       - Integer.  Set rng seed. Empty (default) means don't touch the
%                                     seed.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and coneMosaic) that define the neural
%               compute pipeline for this computation.  This can be useful
%               for a user interested in knowing what needs to be supplied
%               to this.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%                .neuralResponses : dictionary of responses indexed with
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%                .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%                .neuralPipeline  : a struct containing the optics and cone mosaic
%                                   employed in the computation (only returned if
%                                   the parent @neuralResponseEngine object has
%                                   an empty neuralPipeline property)
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags')
%       and are arranged in a matrix of:
%           [instancesNum x mCones x tTimeBins]
%
% See Also:
%     t_neuralResponseCompute

% History:
%    09/26/2020  npc  Wrote it.
%    10/05/2020  dhb  Apply ieParamFormat to varargin for all keys.
%    10/05/2020  dhb  Rename. Work on comments.
%    10/05/2020  dhb  Rewrite to use 'rngSeed' key/value pair.
%    10/17/2020  dhb  Use randomly chosen seed for mosaic compute operation
%                     if rngSeed is set to [].  Save/restore rng state when
%                     an explicit seed is passed.
%                dhb  Just return one response instance in the no noise
%                     case.
%    10/19/2020  dhb  Fix comment to reflect fact that we now return
%                     instancesNum instances in noise free case.

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nrePhotopigmentExcitationsWithNoEyeMovements()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements);

    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngineOBJ = sceneEngine(@sceUniformFieldTemporalModulation);
    testContrast = 0.1;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);
    
    % Compute 16 response instances for a number of different noise flags
    instancesNum = 16;
    noiseFlags = {'random', 'none'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngineOBJ.compute(...
            theTestSceneSequence, ...
            theTestSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );

    % Retrieve the different computed responses
    noiseFreeResponses = theResponses('none');
    randomNoiseResponseInstances = theResponses('random');
%}

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.addParameter('noiseFlags', {'random'});
p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

% Retrieve the response noiseFlag labels and validate them.
noiseFlags = p.Results.noiseFlags;
rngSeed = p.Results.rngSeed;
neuralEngineOBJ.validateNoiseFlags(noiseFlags);

% For each noise flag we generate a corresponing neural response, and all
% neural responses are stored in a dictionary indexed by the noiseFlag label.
% Setup theNeuralResponses dictionary, loading empty responses for now
theNeuralResponses = containers.Map();
for idx = 1:length(noiseFlags)
    theNeuralResponses(noiseFlags{idx}) = [];
end

if (isempty(neuralEngineOBJ.neuralPipeline))
    % Generate the optics
    theOptics = oiCreate(neuralResponseParamsStruct.opticsParams.type, neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    % Generate the cone mosaic
    theConeMosaic = coneMosaicHex(neuralResponseParamsStruct.coneMosaicParams.upsampleFactor, ...
        'fovDegs', neuralResponseParamsStruct.coneMosaicParams.fovDegs, ...
        'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds, ...
        'eccBasedConeDensity', neuralResponseParamsStruct.coneMosaicParams.eccBased, ...
        'eccBasedConeQuantal', neuralResponseParamsStruct.coneMosaicParams.eccBased, ...
        'maxGridAdjustmentIterations', 50 ...
        );
    returnTheNeuralPipeline = true;
else
    % Load the optics from the previously computed neural pipeline
    theOptics = neuralEngineOBJ.neuralPipeline.optics;
    % Load the cone mosaic from the previously computed neural pipeline
    theConeMosaic = neuralEngineOBJ.neuralPipeline.coneMosaic;
    returnTheNeuralPipeline =  false;
end

% Compute the sequence of optical images corresponding to the sequence of scenes
framesNum = numel(sceneSequence);
theListOfOpticalImages = cell(1, framesNum);
for frame = 1:framesNum
    theListOfOpticalImages{frame} = oiCompute(sceneSequence{frame}, theOptics);
end

% Generate an @oiSequence object containing the list of computed optical images
theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);

% Zero eye movements
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
emPaths = zeros(instancesNum, eyeMovementsNum, 2);

% Set rng seed if one was passed. Not clear we need to do this because
% all the randomness is in the @coneMosaic compute object, but it
% doesn't hurt to do so, if we ever choose a random number at this
% level.
if (~isempty(rngSeed))
    oldSeed = rng(rngSeed);
end

% Compute responses for each type of noise flag requested
for idx = 1:length(noiseFlags)
    if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
        % Compute the noise-free response
        % To do so, first save the current mosaic noiseFlag
        lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
        
        % Set the coneMosaic.noiseFlag to 'none';
        theConeMosaic.noiseFlag = 'none';
        
        % Compute noise-free response instances
        theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...   % the emPaths
            'currentFlag', false ...  % no photocurrent
            );
        
        % Restore the original noise flag
        theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
        
    elseif (~isempty(rngSeed))
        % Compute noisy response instances with a specified random noise seed for repeatability
        theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...    % the emPaths
            'currentFlag', false, ...  % no photocurrent
            'seed', rngSeed ...        % random seed
            );
        
    elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
        % Because computeForOISequence freezes noise, if we want
        % unfrozen noise (which is the case if we are here),
        % we have to pass it a randomly chosen seed.
        useSeed = randi(32000,1,1);
        
        % Compute noisy response instances
        theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...   % the emPaths
            'currentFlag', false, ...  % no photocurrent
            'seed', useSeed ...        % random seed
            );
    end
end

% Restore rng seed if we set it
if (~isempty(rngSeed))
    rng(oldSeed);
end

% Temporal support for the neural response
temporalSupportSeconds = theConeMosaic.timeAxis;

% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);
if (returnTheNeuralPipeline)
    dataOut.neuralPipeline.optics = theOptics;
    dataOut.neuralPipeline.coneMosaic = theConeMosaic;
end
end

function p = generateDefaultParams()
% Default params for this compute function
p = struct(...
    'opticsParams', struct(...
    'type', 'wvf human', ...
    'pupilDiameterMM', 3.0 ...
    ), ...
    'coneMosaicParams', struct(...
    'upsampleFactor', 5, ...
    'fovDegs', 0.3, ...
    'timeIntegrationSeconds', 0.1, ...
    'eccBased', false ...
    ) ...
    );
end
