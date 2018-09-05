classdef RCMEoperators < handle
    %RCMEOPERATORSFACTORY RC operators
    
    properties
        RCMEbundle
    end
    
    methods
        
        function this = RCMEoperators(RCMEbundle)
            this.RCMEbundle = RCMEbundle;
        end
        
        function [ Q, H0, Xim, Xid ] = getRCMEoperators( this )
            
            Q = this.operatorQ();
            
            %
            %   Operator H0 in system+RC base.
            %
            H0 = this.operatorH0;
            
            %
            %   Operators Xim and Xid system+RC base.
            %
            [ Xim, Xid ] = this.operatorsXimXid(H0, Q);
        end
        
        function [ Xim, Xid ] = operatorsXimXid( this, H0, Q )
            %OPERATOR Operators Xim and Xid. Cell arrays.
            bundle = this.RCMEbundle;
            %   Temporary go to H0 base
            %
            [eigenvectorMatrix, energyMatrix] = eig(H0);
            Q_eig = cellfun(@(y) eigenvectorMatrix'*y*eigenvectorMatrix, Q, 'UniformOutput', 0);
            H0_eig = energyMatrix;
            %
            %   Create Xim and Xid operators
            %
            Xim = cell(1,bundle.getN());
            Xid = cell(1,bundle.getN());
            for i=1:bundle.getN()
                [Xim{i}, Xid{i}] = this.XimXidTimeIndependant( H0_eig, Q_eig{i} );
            end
            % Go back to initial base
            Xim = cellfun(@(y) eigenvectorMatrix*y*eigenvectorMatrix', Xim, 'UniformOutput', 0);
            Xid = cellfun(@(y) eigenvectorMatrix*y*eigenvectorMatrix', Xid, 'UniformOutput', 0);
        end
        
        
        function [Xim, Xid] = XimXidTimeIndependant( this, H0_eig, Q_eig )
            % One couple of time independant Xim and Xid operators
            bundle = this.RCMEbundle;
            
            Xim = zeros( bundle.getDimAll() );
            Xid = zeros( bundle.getDimAll() );
            for row = 1:bundle.getDimAll()
                for column = 1:bundle.getDimAll()
                    epsilon = H0_eig(row,row) - H0_eig(column,column);
                    % Taylor series
                    if epsilon == 0
                        Xim(row, column) = (1/2)*bundle.spectralDivOmega(epsilon)*(2/bundle.getBeta())*Q_eig(row,column);
                    else
                        Xim(row, column) = (1/2)*bundle.spectral(epsilon)*coth(bundle.getBeta()*epsilon/2)*Q_eig(row,column);
                    end
                    Xid(row, column) = (1/2)*bundle.spectral(epsilon)*Q_eig(row,column);
                end
            end
        end
        
        function H0 = operatorH0(this, varargin)
            %OPERATORH0 Operator H0. 1 matrix matrica.
            
            %System
            H0systemNarys = this.getH0systemPart();
            
            %System-RCs interraction
            H0interactionPart = this.getH0interactionPart();
            
            %RCs
            H0rcPart = this.getH0rcPart();
            
            %All
            H0 = H0systemNarys - H0interactionPart + H0rcPart;
        end
        
        function H0systemPart = getH0systemPart(this)
            bundle = this.RCMEbundle;
            H0s = bundle.hamiltonianSystem();
            H0systemPart = RCMEmath.kron_AIBI(H0s, bundle.identityRC(), 2, bundle.getN()+1);
        end
        
        function H0systemPart = getH0systemPartNoGround(this, withGround)
            bundle = this.RCMEbundle;
            H0s = bundle.hamiltonianSystemNoGround();
            H0systemPart = RCMEmath.kron_AIBI(H0s, bundle.identityRC(), 2, bundle.getN()+1);
        end
        
        
        function [H0interactionPart] = getH0interactionPart(this)
            bundle = this.RCMEbundle;
            
            H0interactionPart = zeros(bundle.getDimAll());
            for RCidx=1:bundle.getN()
                H0interactionPart = H0interactionPart + ...
                    bundle.getEta.*RCMEmath.kron_AIBI(bundle.getOperatoriai{RCidx},...
                    bundle.coordRC(), RCidx+1, bundle.getN()+1);
            end
        end
        
        function H0rcPart = getH0rcPart(this)
            bundle = this.RCMEbundle;
            
            H0rcPart = zeros(bundle.getDimAll());
            for RCidx=1:bundle.getN()
                H0rcPart = H0rcPart + RCMEmath.kron_AIBI(...
                    bundle.identitySystem(), bundle.hamiltonianRC(), RCidx+1, bundle.getN()+1);
            end
        end
        
        function RCsEq = getRCsBoltzmannEquil(this)
            bundle = this.RCMEbundle;
            RCeq = zeros( bundle.getDimRC(), bundle.getDimRC() );
            for i = 1:bundle.getDimRC()
                RCeq( i, i ) = exp( -(bundle.getBeta())*(bundle.getReactDaznis())*(i-1) );
            end
            Z = sum(sum(RCeq));
            RCeq = RCeq./Z;
            RCsEq = RCMEmath.kron_BBB( RCeq, bundle.getN() );
        end
        
        function H0eq = getH0BoltzmannEquil(this, hasGround)
            bundle = this.RCMEbundle;
            H0 = this.operatorH0();
            if hasGround
                H0 = H0(bundle.getDimRCs()+1:end,bundle.getDimRCs()+1:end);
            end
            H0eqState = expm(-bundle.getBeta().*H0);
            H0eqState = H0eqState./trace(H0eqState);
            H0eq = zeros(bundle.getDimAll());
            if hasGround
                H0eq(bundle.getDimRCs()+1:end, bundle.getDimRCs()+1:end) = H0eqState;
            else
                H0eq =H0eqState;
            end
        end
        
        function [ q ] = operatorQ( this )
            %OPERATORA Operator Q (coordinate) for every RC (cell array)
            bundle = this.RCMEbundle;
            q = cell(1,bundle.getN());
            for RCidx=1:bundle.getN()
                q{RCidx} = RCMEmath.kron_AIBI(bundle.identitySystem(), ...
                    bundle.coordRC(), RCidx+1, bundle.getN()+1);
            end
        end
        
        function systemFull = getSystemOperatorFullBasis(this, system)
            bundle = this.RCMEbundle;
            systemFull = RCMEmath.kron_ABBB(system, bundle.identityRC(), bundle.getN()+1);
        end
        
        function fullOperator = combineSystemRCsOperators(this, systemOperator, RCsOperator)
            fullOperator = kron(systemOperator, RCsOperator);
        end
        
    end
    
    methods (Static)
        
        function [ propogator ] = operatorPropogatorOld( Q, H0, Xim, Xid )
            %OPERATORPROPOGATOR Superoperator propogator. Liouville space.
            %Old slow algorithm. Calculates matrix elements.
            
            % Dimension Liouville space
            dimPilnas = size(H0, 1);
            dimLiouville = dimPilnas^2;
            
            % Liuouville and Hilbert dimension corresponance
            % Pvz 1 -> 11, 2 -> 21, 3 -> 12, 4 - >22
            hilbertLiouvilleMap = RCMEutil.createHilbertLiouvilleMap(dimPilnas);
            
            % Helpers
            AXid = cellfun(@(x,y) x*y, Q, Xid, 'UniformOutput',false);
            XidA = cellfun(@(x,y) x*y, Xid, Q, 'UniformOutput',false);
            AXim = cellfun(@(x,y) x*y, Q, Xim, 'UniformOutput',false);
            XimA = cellfun(@(x,y) x*y, Xim, Q, 'UniformOutput',false);
            
            % Creates 
            propogator = spalloc(dimLiouville, dimLiouville, dimLiouville*dimLiouville);
            for row = 1:dimLiouville
                for column = 1:dimLiouville
                    % Indices
                    m = hilbertLiouvilleMap{row}(1);
                    n = hilbertLiouvilleMap{row}(2);
                    j = hilbertLiouvilleMap{column}(1);
                    k = hilbertLiouvilleMap{column}(2);
                    XidIndelis = 0;
                    XimIndelis = 0;
                    for i=1:length(Q)
                        %From Xid
                        XidIndelis = XidIndelis + Q{i}(m,j)*Xid{i}(k,n) - Q{i}(k,n)*Xid{i}(m,j);
                        %From Xim
                        XimIndelis = XimIndelis + Q{i}(m,j)*Xim{i}(k,n) + Q{i}(k,n)*Xim{i}(m,j);
                        
                        if k == n
                            XidIndelis = XidIndelis + AXid{i}(m,j);
                            XimIndelis = XimIndelis - AXim{i}(m,j);
                        end
                        if j == m
                            XidIndelis = XidIndelis - XidA{i}(k,n);
                            XimIndelis = XimIndelis - XimA{i}(k,n);
                        end
                    end
                    %From H0
                    H0indelis = 0;
                    if k == n
                        H0indelis = H0indelis - 1i*H0(m,j);
                    end
                    if j == m
                        H0indelis = H0indelis + 1i*H0(k,n);
                    end
                    % Final
                    propogator(row, column) = XidIndelis + XimIndelis + H0indelis;
                end
            end
        end
        
        function [ propogator ] = operatorPropogator( Q, H0, Xim, Xid )
            %OPERATORPROPOGATOR Superoperator propogator. Liouville space.
            
            % Dimension Liouville space
            dimPilnas = size(H0, 1);
            dimLiouville = dimPilnas^2;
            % Fake density operator, one element will be 1, other 0.
            ro = zeros(dimPilnas, dimPilnas);
            
            propogator = zeros(dimLiouville, dimLiouville);
            
            % For every column
            for propColIndex = 1:dimLiouville
                % Marks column
                ro(propColIndex) = 1;
                
                % Calculates column
                propCol = -1i*RCMEmath.commutator(H0, ro);
                for i=1:length(Q)
                    propCol = propCol - RCMEmath.commutator(Q{i}, RCMEmath.commutator(Xim{i}, ro)) + ...
                        RCMEmath.commutator(Q{i}, RCMEmath.antiCommutator(Xid{i}, ro));
                end
                % Saves
                propCol = propCol(:);
                propogator(:, propColIndex) = propCol;
                
                ro(propColIndex) = 0;
            end
            propogator = sparse(propogator);
        end
        
    end
    
end

