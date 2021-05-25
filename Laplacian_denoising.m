function [Q_TV] = Laplacian_denoising(Q_TV,direction,numb_filtered)

% input:
% Q_TV: structure containing original Laplacian (L), original eigenvectors (V,W) and eigenvalues (D)
% direction: string containing the type of filtering 
%            lp: only smaller eigenvalues and associated eigenvectors; 
%            hp: only smalllarger eigenvalues and associated eigenvectors; 
%            lphp: our proposed Laplacian denoising 
% numb_filtered: number of components (eigenvalues and eigenvectors which are kept to built the new Laplacian)

% output
% Q_TV: structure containing original matrices and  Laplacian after denoising



for t=1:Q_TV.Ntrials % cycle on trials
    for w=1:Q_TV.T   % cycle on time windows time window
        
                switch direction
                    case 'lp'% filtering with sameller eigenvectors
                        Q_TV.D_filt(:,:,w,t)=Q_TV.D(1:numb_filtered,1:numb_filtered,w,t);
                        Q_TV.V_filt(:,:,w,t)=Q_TV.V(:,1:numb_filtered,w,t);
                        Q_TV.W_filt(:,:,w,t)=Q_TV.W(:,1:numb_filtered,w,t);

                    case 'hp' % filtering with larger eigenvectors
                        Q_TV.D_filt(:,:,w,t)=Q_TV.D(end-numb_filtered+1:end,end-numb_filtered+1:end,w,t);
                        Q_TV.V_filt(:,:,w,t)=Q_TV.V(:,end-numb_filtered+1:end,w,t);
                        Q_TV.W_filt(:,:,w,t)=Q_TV.W(:,end-numb_filtered+1:end,w,t);
    
                    case 'lphp' % proposed denoising algorithm (smaller+larger)                      
                        Q_TV.D_filt(:,:,w,t)=Q_TV.D([1:numb_filtered+1, end-numb_filtered:end],[1:numb_filtered+1, end-numb_filtered:end],w,t);
                        Q_TV.V_filt(:,:,w,t)=Q_TV.V(:,[1:numb_filtered+1, end-numb_filtered:end],w,t);
                        Q_TV.W_filt(:,:,w,t)=Q_TV.W(:,[1:numb_filtered+1, end-numb_filtered:end],w,t);

                        
                end
                Q_TV.L_filt(:,:,w,t)= Q_TV.W_filt(:,:,w,t)*Q_TV.D_filt(:,:,w,t)*Q_TV.V_filt(:,:,w,t)';
                %Q_TV.A_filt(:,:,w,t)=-(Q_TV.L_filt(:,:,w,t)-Q_TV.L_filt(:,:,w,t).*diag(ones(Q_TV.N,1)));
                    
        
    end
end
