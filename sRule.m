classdef sRule < handle
    %SRULE Summary of this class goes here
    %   Detailed explanation goes here
    
    %properties ( Access = private )
    properties
        ruleBook = struct(); % a struct containing struct containing function handle
        %           ex: struct(,...
        %                   'opDFT',struct(...
        %                               'unit',@,...
        %                               'delta',@),...
        %                   'opBlah',struct(...
        %                               'label',@,...
        %                               'unit',@)...
        %                   )
        %
    end
    
    methods
        function r = sRule(varargin)
            % no argument
            % or one argument of string containing the class name
            if nargin == 1
                opClassName = varargin{1};
                recognizedFields = sRule.recognizedContainerFields();
                
                for ind = 1:length(recognizedFields)
                    if ismethod(opClassName,strcat('default_',recognizedFields{ind},'_handler')) % if default handler is defined in the op
                        r.setRule(recognizedFields{ind},eval(strcat(opClassName,'.default_',recognizedFields{ind},'_handler')));
                    end
                end
            elseif nargin > 1
                error('no or 1 input for sRule constructor')
            end
        end
        
%         function res = getRuleList(sRuleObject,useDefault)
%             % return a struct containing function handle on each field.
%             % operatorName could be a string for the operator name 
%             % or a cell array of string of operator name
%             
%             % useDefault is set to false
%             if nargin < 3, useDefault = false; end
%             
%             % if nothing specified, return the entire rule book
%             if nargin == 1
%                 res = sRuleObject.ruleBook;
%                 return;
%             end
%             
%             % wrap operator name in cell
%             if ischar(operatorName)
%                 operatorName = {operatorName};
%             end
%             
%             %initialize rule book
%             res = struct();
%             
%             if useDefault == 0
%                 fn = fieldnames(sRuleObject.ruleBook); %list of operator names in rule book
%                 while (~isempty(operatorName)) % while the to-process list is not empty
%                     thisOperatorHandled = false; % the operator is not handled (yet)
%                     
%                     for ind = 1:length(fn) % for every operator rule in the rule book
%                         if strcmp(operatorName{1},fn{ind}) % if the the operator is found on the rule book
%                             res.(fn{ind}) = sRuleObject.ruleBook.(operatorName{1}); % add the rule
%                             thisOperatorHandled = true; % mark the operator has handled
%                             break;
%                         end
%                     end
%                     if ~thisOperatorHandled % if not handled
%                         fn = sRule.recognizedContainerFields;
%                         for ind = 1:length(fn)
%                             %get the default handler of the fields from the operator class
%                             res.(operatorName{1}).(fn{ind}) = eval(strcat(operatorName{1},'.default_',fn{ind},'_handler'));
%                         end
%                         thisOperatorHandled = true;
%                     end
%                     operatorName(1) = [];
%                 end
%                 
%             elseif useDefault == 1
%                 fn = sRule.recognizedContainerFields;
%                 for ind = 1:length(fn)
%                     res.(operatorName).(fn{ind}) = eval(strcat(operatorName,'.default_',fn{ind},'_handler'));
%                 end
%             else
%                 error('invalid useDefault argument')
%             end
%         end

        function res = setRule(sRuleObject,varargin)
            for ind = 1:2:length(varargin)
                field = varargin{1}; fHandle = varargin{2};
                if any(ismember(fieldnames(sRuleObject.ruleBook),field))
                    sRuleObject.ruleBook = rmfield(sRuleObject.ruleBook,field);
                end
                sRuleObject.ruleBook.(field) = fHandle;
                varargin(1) = []; varargin(1) = [];
            end
            res = 1;
        end
    end
    
    methods ( Static )
%         function res = recognizedOperators
%             % returns a cell array of strings
%             res = {...
%                 'opDFT',...
%                 'opCurvelet',...
%                 'opWavelet'... % TODO: add more operators
%             };
%         end
%         
%         function res = isRecognizedContainerField(theFieldName)
%             res = false;
%             fn = sRule.recognizedContainerFields;
%             for ind = 1:length(fn)
%                 if strcmpi(theFieldName,fn{ind})
%                     res = true;
%                     return;
%                 end
%             end
%         end
        function res = recognizedContainerFields
           % returns a cell array of strings
            res = {...
                'full',...
                'varName',...
                'varUnits',...
                'origin',...
                'delta',...
                'label',...
                'size',...
                'unit'...
            };
        end
        
        function fh = default_handler(fn)
            fh = eval(strcat('@(header,op,mode)header.',fn));
        end
    end
end

