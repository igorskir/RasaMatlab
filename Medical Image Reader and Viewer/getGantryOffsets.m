function [gantryOffset1, gantryOffset2, varargout] = ...
    getGantryOffsets(gantryModel, patientOrientation, IOP, voxDim1, voxDim2)

%this function is used with interfile formats

%scanner-dependent corrections are different between mCT and mMR
if strcmpi(gantryModel(end-2:end),'mCT') || strcmpi(gantryModel,'1104')
    
    %empirical corrections for mCT gantry offsets use patientOrientation
    %NOTE: IOP may not be accurate at this point for mCT (feet first cases)
    %   and will be changed in readImages function later
    
    %gantryOffset1 = -1 * (0.107 + voxDim1/2);
    %gantryOffset2 = -1 * (-0.417 + voxDim2/2);
    %gantryOffset1 = -1 * (0.10602 + voxDim1/2);
    %gantryOffset2 = -1 * (-0.41599 + voxDim2/2);
    switch patientOrientation
        %define parameters for IPP for each case
        case 'HFS'
            if all(IOP==[1;0;0;0;1;0])
                %IOP does not need to be changed
                gantryOffset2 = -1 * ...
                    (-0.416966604863952 + 0.500486395175824 * voxDim2);
            else
                error('Check this')
            end
        case 'FFS'
            if all(IOP==[1;0;0;0;1;0])
                IOP = [-1; 0; 0; 0; 1; 0];
                gantryOffset2 = -1 * ...
                    (-0.416966604863952 + 0.500486395175824 * voxDim2);
            else
                error('Check this')
            end
        case 'HFP'
            if all(IOP==[-1;0;0;0;-1;0])
                %IOP does not need to be changed
                gantryOffset2 = -1 * ...
                    (-0.416966604863952 + 0.500486395175824 * voxDim2);
            else
                error('Check this')
            end
        case 'FFP'
            if all(IOP==[-1;0;0;0;-1;0])
                IOP = [1; 0; 0; 0; -1; 0];
                gantryOffset2 = -1 * ...
                    (-0.416966604863952 + 0.500486395175824 * voxDim2);
            else
                error('Check this')
            end
        case 'HFDL'
            if all(IOP==[0;-1;0;1;0;0])
                %IOP does not need to be changed
                gantryOffset2 = -1 * ...
                    (-0.417058105013725 + 0.500501372513394 * voxDim2);
            else
                error('Check this')
            end
        case 'FFDL'
            if all(IOP==[0;-1;0;1;0;0])
                IOP = [0; 1; 0; 1; 0; 0];
                gantryOffset2 = -1 * ...
                    (-0.417058105013725 + 0.500501372513394 * voxDim2);
            else
                error('Check this')
            end
        case 'HFDR'
            if all(IOP==[0;1;0;-1;0;0])
                %IOP does not need to be changed
                gantryOffset2 = -1 * ...
                    (-0.416966604863952 + 0.500486395175824 * voxDim2);
            else
                error('Check this')
            end
        case 'FFDR'
            if all(IOP==[0;1;0;-1;0;0])
                IOP = [0; -1; 0; -1; 0; 0];
                gantryOffset2 = -1 * ...
                    (-0.417058105013725 + 0.500501372513394 * voxDim2);
            else
                error('Check this')
            end
        otherwise
            error('Patient orientation is not recognized')
    end
    temp = cross(IOP(1:3),IOP(4:6));
    gantryOffset1 = -1 * ((0.106486098507667 - temp(3) * 0.001449503371619) + ...
        (0.500149233209423 + temp(3) * 0.0003371619664019332) * voxDim1);
elseif strcmpi(gantryModel(end-2:end),'mMR') || strcmpi(gantryModel,'2008')
    
    %empirical corrections for mMR gantry offsets depend on IOP
    
    if abs(IOP(1))==1 || abs(IOP(2))==1
        gantryOffset1 = (-1 * (1.02334 + voxDim2/2) * IOP(4)) + ... % HFDL, HFRL, FFDL, or FFDR
            (-1 * (-0.1450 + (voxDim1/2 * IOP(5))) * IOP(1) * IOP(5)); % HFS, FFS, HFP, or FFP
        gantryOffset2 = (-1 * (-0.1450 + voxDim1/2) * IOP(2)) + ... % HFDL, HFRL, FFDL, or FFDR
            (-1 * (1.02334 + (voxDim2/2 * IOP(5))) * IOP(5) * IOP(5)); % HFS, FFS, HFP, or FFP
        %gantryOffset3 = -1 * 0.1724;
        temp = cross(IOP(1:3),IOP(4:6));
        gantryOffset3 = -1 * temp(3) * (0.00025 + temp(3) * 0.1719);
        varargout{1} = gantryOffset3; %mCT does not require gantryOffset3 output
    else
        error('Patient orientation is not recognized')
    end
end

end