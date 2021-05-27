% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           This function adds one additional incoherent layer 
%                         and a coherent layer stack calculated by TMM beforehand
% #INPUT:                 Rf1: Reflection of initial architecture in forward direction (tensor)
%                         Af1: Absorption of initial architecture in forward direction (tensor)
%                         Tf1: Transmission of initial architecture in forward direction (tensor)
%                         Rb1: Reflection of initial architecture in backward direction (tensor)
%                         Ab1: Absorption of initial architecture in backward direction (tensor)
%                         Tb1: Transmission of initial architecture in backward direction (tensor)
%                         Aincoherent: Absorption of the added incoherent layer (array)
%                         Rf2: Reflection of effective interface after the incoherent layer in forward direction (tensor)
%                         Af2: Absorption of effective interface after the incoherent layer in forward direction (tensor)
%                         Tf2: Transmission of effective interface after the incoherent layer in forward direction (tensor)
%                         Rb2: Reflection of effective interface after the incoherent layer in forward direction (tensor)
%                         Ab2: Absorption of effective interface after the incoherent layer in forward direction (tensor)
%                         Tb2: Transmission effective interface after the incoherent layer in forward direction (tensor)
% #OUTPUT:                Rf: Reflection of enhanced architecture in forward direction (array)
%                         Af: Absorption of enhanced architecture in forward direction (array)
%                         Tf: Transmission of enhanced architecture in forward direction (array)
%                         Rb: Reflection of enhanced architecture in backward direction (array)
%                         Ab: Absorption of enhanced architecture in backward direction (array)
%                         Tb: Transmission of enhanced architecture in backward direction (array)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          
%                         
% -----------------------------------------------------------

function [Rf,Af,Tf,Rb,Ab,Tb] = AddIncoherentCoherent(Rf1,Af1,Tf1,Rb1,Ab1,Tb1,Aincoherent,Rf2,Af2,Tf2,Rb2,Ab2,Tb2)

    Rf = zeros(size(Rf1));
    Af = zeros(size(Af1,1)+size(Af2,1)+1,size(Af1,2),size(Af1,3));
    Tf = zeros(size(Tf1));
    Rb = zeros(size(Rf1));
    Ab = zeros(size(Af1,1)+size(Af2,1)+1,size(Af1,2),size(Af1,3));
    Tb = zeros(size(Tf1));
    for idx=1:size(Rf1,3)
        RefIncoherentForward = diag(1-Aincoherent(:,idx))*Rf2(:,:,idx)*diag(1-Aincoherent(:,idx));
        RefIncoherentBackward = diag(1-Aincoherent(:,idx))*Rb1(:,:,idx)*diag(1-Aincoherent(:,idx));
        
        % Reflection
        Rf(:,:,idx) = Rf1(:,:,idx)+... % First reflection at initial interface forward
            Tf1(:,:,idx)*RefIncoherentForward*Tb1(:,:,idx); % First reflection at interface to the next layer(s) forward
        Rb(:,:,idx) = Rb2(:,:,idx)+... % First reflection at interface to the next layer(s) backward
            Tb2(:,:,idx)*RefIncoherentBackward*Tf2(:,:,idx); % First reflection at initial interface backward
        for j=1:1000 % Reflection with added layer forward
            x = Tf1(:,:,idx)*(RefIncoherentForward*Rb1(:,:,idx))^j*RefIncoherentForward*Tb1(:,:,idx); % Multiple reflections between initial interace and interface to the next layer(s)
            Rf(:,:,idx) = Rf(:,:,idx) + x;
            if max(x(:))<0.0001
                break;
            end
        end
        for j=1:1000 % Reflection with added layer backward
            x = Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx))^j*RefIncoherentBackward*Tf2(:,:,idx); % Multiple reflections between initial interace and interface to the next layer(s)
            Rb(:,:,idx) = Rb(:,:,idx) + x;
            if max(x(:))<0.0001
                break;
            end
        end
        
        % Transmission
        Tf(:,:,idx) = Tf1(:,:,idx)*diag(1-Aincoherent(:,idx))*Tf2(:,:,idx);
        Tb(:,:,idx) = Tb2(:,:,idx)*diag(1-Aincoherent(:,idx))*Tb1(:,:,idx);
        for l=1:1000 % Transmission with added layer forward
            z = Tf1(:,:,idx)*(diag(1-Aincoherent(:,idx))*Rf2(:,:,idx)*diag(1-Aincoherent(:,idx))*Rb1(:,:,idx))^l*diag(1-Aincoherent(:,idx))*Tf2(:,:,idx);
            Tf(:,:,idx) = Tf(:,:,idx) + z;
            if max(z(:))<0.0001
                break;
            end
        end
        for l=1:1000 % Transmission with added layer backward
            z = Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx))^l*diag(1-Aincoherent(:,idx))*Tb1(:,:,idx);
            Tb(:,:,idx) = Tb(:,:,idx) + z;
            if max(z(:))<0.0001
                break;
            end
        end
        
        % Absorption
        Af(1:size(Af1,1),:,idx) = Af1(:,:,idx)+...
            (Tf1(:,:,idx)*RefIncoherentForward*Ab1(:,:,idx)')'; % Absorption in initial layers forward
        Af(size(Af1,1)+1,:,idx) = Tf1(:,:,idx)*Aincoherent(:,idx)+...
            Tf1(:,:,idx)*diag(1-Aincoherent(:,idx))*Rf2(:,:,idx)*Aincoherent(:,idx); % Absorption in added incoherent layer forward
        Af(size(Af1,1)+2:end,:,idx) = (Tf1(:,:,idx)*diag(1-Aincoherent(:,idx))*Af2(:,:,idx)')'; % Absorption in added coherent layer stack forward
        
        Ab(1:size(Af1,1),:,idx) = (Tb2(:,:,idx)*diag(1-Aincoherent(:,idx))*Ab1(:,:,idx)')'; % Absorption in initial layers backward
        Ab(size(Af1,1)+1,:,idx) = Tb2(:,:,idx)*Aincoherent(:,idx)+... % Absorption in added incoherent layer backward
            Tb2(:,:,idx)*diag(1-Aincoherent(:,idx))*Rb1(:,:,idx)*Aincoherent(:,idx);
        Ab(size(Af1,1)+2:end,:,idx) = Ab2(:,:,idx)+...
            (Tb2(:,:,idx)*RefIncoherentBackward*Af2(:,:,idx)')'; % Absorption in added coherent layer stack backward
        
        for j=1:1000 % Absorption forward
            x = Tf1(:,:,idx)*(RefIncoherentForward*Rb1(:,:,idx))^j*RefIncoherentForward*Ab1(:,:,idx)';
            y = Tf1(:,:,idx)*(RefIncoherentForward*Rb1(:,:,idx))^j*Aincoherent(:,idx)+...
                Tf1(:,:,idx)*(RefIncoherentForward*Rb1(:,:,idx))^j*diag(1-Aincoherent(:,idx))*Rf2(:,:,idx)*Aincoherent(:,idx);
            z = Tf1(:,:,idx)*(RefIncoherentForward*Rb1(:,:,idx))^j*diag(1-Aincoherent(:,idx))*Af2(:,:,idx)';
            Af(1:size(Af1,1),:,idx) = Af(1:size(Af1,1),:,idx) + x';
            Af(size(Af1,1)+1,:,idx) = Af(size(Af1,1)+1,:,idx) + y';
            Af(size(Af1,1)+2:end,:,idx) = Af(size(Af1,1)+2:end,:,idx) + z';
            if max(x(:))<0.0001 && max(y(:))<0.0001 && max(z(:))<0.0001
                break;
            end
        end
        for j=1:1000 % Absorption backward
            x = Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx)*diag(1-Aincoherent(:,idx)))^j*Ab1(:,:,idx)';
            y = Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx))^j*Aincoherent(:,idx)+...
                Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx))^j*diag(1-Aincoherent(:,idx))*Rb1(:,:,idx)*Aincoherent(:,idx);
            z = Tb2(:,:,idx)*(RefIncoherentBackward*Rf2(:,:,idx))^j*RefIncoherentBackward*Af2(:,:,idx)';
            Ab(1:size(Ab1,1),:,idx) = Ab(1:size(Ab1,1),:,idx) + x';
            Ab(size(Ab1,1)+1,:,idx) = Ab(size(Ab1,1)+1,:,idx) + y';
            Ab(size(Ab1,1)+2:end,:,idx) = Ab(size(Ab1,1)+2:end,:,idx) + z';
            if max(x(:))<0.0001 && max(y(:))<0.0001 && max(z(:))<0.0001
                break;
            end
        end
    end
end