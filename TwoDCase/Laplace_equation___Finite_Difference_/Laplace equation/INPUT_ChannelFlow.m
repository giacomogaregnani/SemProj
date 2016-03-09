function [channel_Dimension,obstruction_Dimension,grid_values,InletVel] = INPUT_ChannelFlow()

disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('%                    %')
disp('%                    %')
disp('%%%%%%%%%%%%%%%%%%%%%%')
ChW = input('Channel Width >> ');ChL = input('Channel Length >> ');ISz = input('Inlet Size >> ');
channel_Dimension = [ChL ; ChW ; ISz];
disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('Please input :-')
disp('    1 for channel flow across cylinder.')
disp('    2 for channel flow across rectangle.')
disp('    3 for unobstructed channel flow.')
for count = 1:1000
    value = input('Your input please (def. = 1) >> ');
    if isempty(value)
        value   = 1; % This is the value of the default
        default = 1; % User accepts default value
    else
        default = 0; % User wishes to input custom value
    end
    if default   == 1
        break
    else % Then, user input is checked for correctness
        if value == 1 || value == 2 || value == 3
            break
        else
            disp('Please input the correct values "or" use default')
        end
    end
end
Obstruction = value;clear value
disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('%                    %')
disp('%                    %')
disp('%%%%%%%%%%%%%%%%%%%%%%')
if Obstruction == 1
    disp('You chose Cyl. as the obst. ....')
    disp('Enter 1 for circular cyl. and 2 for ellipsoidal cyl.')
    disp('Circ. type is default')
    CylType = input('Your choice of cyl. (def. = 1) >> ');
    if isempty(CylType)
        CylType = 1;
    end
    if CylType == 1
        disp('You chose Circ. Cyl. as the obst. ....')
        CylDia = input('Cyl. dia. >> ');
        if CylDia > (3/4)*ChW
            for count = 1:10000
                fprintf('Cyl. Dia must be ( <= %d)',(3/4)*ChW )
                CylDia = input('Cyl. dia. >> ');
                if CylDia <= (3/4)*ChW
                    break
                end
            end
        end
        CylDist = input('Cyl. Dist. from x = 0 >>');
        if CylDist <= ISz || CylDist >= (ChL - CylDia)
            for count = 1:10000
                fprintf('Cyl. Dist. must be ( > %d ) & ( < %d ) \n',ISz,ChL - CylDia)
                CylDist = input('Cyl. Dist. from x = 0 >>');
                if CylDist > ISz && CylDist < (ChL - CylDia)
                    break
                end
            end
        end
        fprintf('Cyl. Off. must be ( < %d)\n',1.25*(ChW - CylDia)/2)
        CylOff = input('Your value for Cyl. Off. (def. = 0) >> ');
        default = 0;
        if isempty(CylOff)
            CylOff = 0; % Assign default value
            default = 1;            
        end
        if default == 0
            if CylOff >= 1.25*(ChW - CylDia)/2
                for count = 1:10000
                    fprintf('Cyl. Off. must be ( < %d)\n',1.25*(ChW - CylDia)/2)
                    CylOff = input('Your value for Cyl. Off. (def. = 0) >> ');
                    if isempty(CylOff)
                        CylOff = 0;
                    end
                    if CylOff < 1.25*(ChW - CylDia)/2
                        break
                    end
                end
            end
        end
        clear default count
        obstruction_Dimension = [1 ; CylDist ; CylDia ; CylDia ; 0 ; CylOff];
        clear CylDia CylOff
    elseif CylType == 2
        disp('You chose Ellip. Cyl. as the obst. ....')
        CylDiaMajor           = input('Cyl. maj. dia. >> ');
        if CylDiaMajor == 0 || CylDiaMajor > (3/4)*ChW
            for count = 1:10000
                fprintf('Cyl. Mj. Dia must be ( > 0 and <= %d)',(3/4)*ChW )
                CylDiaMajor = input('Cyl. dia. >> ');
                if CylDiaMajor > 0 && CylDiaMajor <= (3/4)*ChW
                    break
                end
            end
        end
        CylDiaMinor = input('Cyl. min. dia. (def.: r_maj = r_min) >> ');
        default = 0;
        if isempty(CylDiaMinor)
            CylDiaMinor = CylDiaMajor;
            default = 1;
        end
        if default == 0
            if CylDiaMinor == 0 || CylDiaMinor < 0
                for count = 1:10000
                    disp('Cyl. Min. Dia cannot be zero and -ve ..')
                    CylDiaMinor = input('Cyl. min. dia. (def.: r_maj = r_min) >> ');
                    if isempty(CylDiaMinor)
                        CylDiaMinor = CylDiaMajor;
                    end
                    if CylDiaMinor <= CylDiaMajor && CylDiaMinor > 0
                        break
                    end
                end
            elseif CylDiaMinor > CylDiaMajor
                for count = 1:10000
                    disp('Min. dia. must be > 0 and < maj. dia.')
                    CylDiaMinor = input('Cyl. min. dia. (def.: r_maj = r_min) >> ');
                    if isempty(CylDiaMinor)
                        CylDiaMinor = CylDiaMajor;
                    end
                    if CylDiaMinor <= CylDiaMajor && CylDiaMinor > 0
                        break
                    end
                end
            end
        end
        clear default
        CylDist = input('Cyl. Dist. from x = 0 >>');
        if CylDist <= (ISz + CylDiaMajor/2 - CylDiaMinor/2) || CylDist >= (ChL - CylDiaMajor)
            for count = 1:10000
                fprintf('Cyl. Dist. must be ( > %2.4f ) & ( < %2.4f ) \n',ISz + CylDiaMajor/2 - CylDiaMinor/2,ChL - CylDiaMajor)
                CylDist = input('Cyl. Dist. from x = 0 >>');
                if CylDist > (ISz + CylDiaMajor/2 - CylDiaMinor/2) && CylDist < (ChL - CylDiaMajor)
                    break
                end
            end
        end
        CylAngle = input('Cyl. inclination about horz., in deg. (def. = 0) >> ');
        default = 0;
        if isempty(CylAngle)
            CylAngle = 0;
            default = 1;
        end
        if default == 0
            if CylAngle > 360
                CylAngle = +mod(CylAngle,360);
            elseif CylAngle < -360
                CylAngle = -mod(sqrt(CylAngle^2),360);
            end
        end
        clear default
        CylOff = input('Cyl. Offset from Ch. axis (def. = 0) >> ');
        default = 0;
        if isempty(CylOff)
            CylOff = 0;
            default = 0;
        end
        if default == 0
            if CylOff >= 0.75*(ChW - CylDiaMajor)/2
                for count = 1:10000
                    fprintf('Cyl. Off. must be ( < %2.4f)\n',0.75*(ChW - CylDiaMajor)/2)
                    CylOff = input('Your value for Cyl. Off. (def. = 0) >> ');
                    if isempty(CylOff)
                        CylOff = 0;
                    end
                    if CylOff < 0.75*(ChW - CylDiaMajor)/2
                        break
                    end
                end
            end
        end
        clear default
        obstruction_Dimension = [2 ; CylDist ; CylDiaMajor ; CylDiaMinor ; CylAngle ; CylOff];
        clear CylDiaMajor CylDiaMinor CylAngle CylOff CylDist
    end
    clear CylType
elseif Obstruction == 2
    disp('You chose Rect. as the obstr. ....')
    RecDist                   = input('Distance from inlet >> ');
    if isempty(RecDist)
        RecDist = channel_Dimension(1,1)/2;
    end
    RecBreadth                = input('Rect. Breadth >> ');
    if isempty(RecBreadth)
        RecBreadth = channel_Dimension(2,1)/2;
    end
    RecLength                 = input('Rect. Length >> ');
    if isempty(RecLength)
        RecLength = RecBreadth;
    end
    RecAngle                  = input('Rect. inclination with horz. in deg. (def. = 0) >> ');
    if isempty(RecAngle)
        RecAngle = 0;
    end
    RecOff                    = input('Rect. Offset from Channel axis (def. = 0) >> ');
    if isempty(RecOff)
        RecOff = 0;
    end
    obstruction_Dimension     = [3 ; RecDist ; RecBreadth ; RecLength ; RecAngle ; RecOff];
    clear RecBreadth RecLength RecAngle RecOff
elseif Obstruction == 3
    obstruction_Dimension     = [4 ; NaN ; NaN ; NaN ; NaN ; NaN];
end
disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('%                    %')
disp('%                    %')
disp('%%%%%%%%%%%%%%%%%%%%%%')
Xgrid = input('Grid Spacing along X-axis >>');
fprintf('Maximum Y grid spacing allowed is %2.4f \n',ISz/2)
Ygrid = input('Grid Spacing along Y-axis >>');
grid_values = [Xgrid Ygrid];
disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('%                    %')
disp('%                    %')
disp('%%%%%%%%%%%%%%%%%%%%%%')
InletVel = input('Please input the inlet velocity >>');
disp('%%%%%%%%%%%%%%%%%%%%%%')
disp('%                    %')
disp('%                    %')
disp('%%%%%%%%%%%%%%%%%%%%%%')