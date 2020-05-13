% spectral_cluster:
% Given the transition matrix for a state space with N states,
% reduce the N states into Mc clusters, 
% using eigenvectors of the transition matrix.
% Briefly, split the N states into 2 clusters using the largest eigenvector
% (with non-unitary eigenvalue), then split one of those two clusters into
% 2 to create 3 clusters, by choosing the cluster with the largest
% variance, and continue until Mc clusters exist. 
%%%%%%INPUTS:
%1. ptrans: NxN matrix of transition probabilities. row i, col j defines
%probability of transitioning from state i to state j. Sum over all columns in a row thus adds to 1.
%Matrix is not symmetric (usually), because transition from i->j does not equal transition from j->i
%2. peq: Nx1 vector of equilibrium probabilities for each state. From an MD
%trajectory, this will be the number of times you are found in that state,
%normalized by number of configurations total.
%3. lag_time: float value (1x1). The amount of time between consequtive
%configurations in the trajectory used to construct the transition matrix.
%4. Mc: integer (1x1). Number of clusters to reduce your N states into.
%%%%OUTPUTS:
%1.Ptransnew: McxMc. The new transition probabitily matrix, based on re-assigning
%your N states in to Mc states.
%2.peqnew: Mcx1. The new equilibrium probability for being in each of your Mc
%states.
%3.clust_reassignments: Mcx2 matrix. First column is each of the 1:N states, Second
%column is the new state it was re-assigned to, numbered 1:Mc. 
%4. eigvalue: Nx1: eigenvalues of the original transition matrix.
%5. evnew: Mcx1: eigenvalues of the new, reassigned states transition
%matrix
%6. taunew: Mc-1x1: Relaxation times controlling the full system dynamics.
%7. K: McxMc: matrix of rates for transitioning between states, units of
%(lag_time)^-1.
%8. time_states: Mcx1: average time spent in each state (units lag_time). 
%%%%%


function[Ptransnew, peqnew, clust_reassignments, eigvalue, evnew, taunew, K, time_states ]=spectral_cluster(ptrans,peq, lag_time, Mc)

%number of current clusters
[Nc, ig]=size(ptrans);

%ENFORCE DETAILED BALANCE (DB)
%peq_ptrans=diag(peq)*ptrans
%ps=0.5*(peq_ptrans+peq_ptrans');
%ptransDB=diag(1./peq)*ps;

%rate matrix
K=logm(ptrans)/lag_time;%because P_trans=exp(K*lag_time), where it is a matrix exponential.


[v,d,w]=eig(ptrans);

[eigvalue, index]=sort(real(diag(d)),'descend');%put eigenvalues in order, save re-ordering indices. 

lambda=log(eigvalue)/lag_time;%eigenvalues of rate matrix
tau=-1./lambda(2:end); %timescales of transitions, skip first e-value as it is the steady state (tau=0).
f1=figure(1)
axes1=axes('Parent',f1,'LineWidth',2,'FontWeight','demi','FontSize',40);
hold(axes1);
plot([2:1:(length(tau)+1)], tau,'ro-','MarkerSize',5)
xlabel('Eigenvalue Index')
ylabel('\tau=-lagtime/log(eig(ptrans))')
legend('Gaps in this spectrum indicate a good number of clusters to choose')



%The first eigenvector corresponding to evalue=1 is the equilibrium
%distribution of states. w(:,1)/sum(w(:,1))=peq (with some numerical
%error).

%use largest eigenvalues (not=1) to split into new clusters.
%to create Mc clusters, need to split from 1:2:3:4:Mc, which is Mc-1
%splits.
nsplit=Mc-1;

states=ones(Nc,1);%assignments of where original clusters (each row of ptrans) go.
ivec=zeros(Nc, Mc);
%m is also the number of states created so far.
for m=1:1:nsplit
    currev=eigvalue(m+1); %current eigenvalue
    currind=index(m+1); %index of current eigenvalue/eigenvector pair
    evec=real(w(:,currind));
    %evaluate the variance of all current states, split the one with
    %largest variance
    vmax=-1;
    for i=1:1:m
        inds=find(states==i);
        v(i)=var(evec(inds));%evaluate variance of these states in the eigenvector
        if(v(i)>vmax)
            vmax=v(i);
            ssplit=i;
        end
    end
    %state to split
    ssplit
    inds=find(states==ssplit);
    svec=evec(inds);%these are eignevector elements only for this state.
    smax=max(svec);
    smin=min(svec);
    sig=(svec-smin)/(smax-smin);%values between 0 and 1
    s1=find(sig>0.5);
    s2=find(sig<0.5);
    %now change the state identities of the elements in s1 and s2, mapped
    %back to ivec, into two states.
    newstate=m+1;
    sold=inds(s1);
    snew=inds(s2);
    states(sold)=ssplit;
    states(snew)=newstate;
    states
end

peqnew=zeros(Mc, 1);%new equilibrium vector
Ptransnew=zeros(Mc, Mc);%new transition matrix
for m=1:1:Mc
    iv=find(states==m);
    peqnew(m)=sum(peq(iv));
end
%update transition matrix
for m=1:1:Mc
    iv=find(states==m);
    
    for i=m:1:Mc
        iv2=find(states==i);
        Ptransnew(m,i)=sum(sum(ptrans(iv, iv2)));
        Ptransnew(i, m)=sum(sum(ptrans(iv2, iv)));
    end
end
%renormalize transition matrix such that for each row, columns sum to 1.
for m=1:1:Mc
    sval=sum(Ptransnew(m,:));%this should just be the number of states in the cluster.
    Ptransnew(m,:)=Ptransnew(m,:)/sval;
end

origstates=[1:1:Nc]';
newstates=states;
clust_reassignments=[origstates, newstates];

[v,d,w]=eig(Ptransnew);
evnew=sort(real(diag(d)),'descend');
lambda=log(evnew)/lag_time;
taunew=-1./lambda(2:end);%relaxation times
K=logm(Ptransnew)/lag_time; %rates for transitions (units of (lag_time)^-1)
kdiag=diag(K);
time_states=-1./kdiag;%average time spent in each state.

    
