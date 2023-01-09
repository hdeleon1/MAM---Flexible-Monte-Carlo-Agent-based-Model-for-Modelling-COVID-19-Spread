%function [Who_ill]=MAM_deleon(N,R,who_vac)
%a2_infec - the chance to get infected after 21 days from the first dose
%a3_infec - the chance to get infected after 180 days from the first dose
N=1.1*1e4;
load('R.csv')
ratio=(1.25./sqrt(R));
ages=[100,82,73,61,48,34,0];
Dages=abs(diff(ages));
load('who_vec3.mat')
load('Fit.mat')


Who_vac=who_vac;

x0=linspace(0,5e3,5e3+1);
y0=linspace(0,5e3,5e3+1);
r0=sqrt(x0.^2+y0.^2);
patient_zero_x=x0(ceil(2e3*rand(1,1)));
patient_zero_y=y0(ceil(2e3*rand(1,1)));
days=[];
vec_days=0;
number_of_people=ceil(rand(1)*5);
n_not_ill=N-1;
n_ill=N-n_not_ill;
M_ill_x=patient_zero_x;
M_ill_y=patient_zero_y;

Who_free=(1:N);
Who_free=Who_free(randperm(length(Who_free)));
Who_free=Who_free(1:0.95*round(length(Who_free)));
Who_rec=[];
Who_ill=[ceil(rand(1,1)*N),0];
Who_not_ill=[Who_free];
Who_not_ill=Who_not_ill(Who_not_ill~=Who_ill(1,1));

n_not_ill=length(Who_not_ill);

M_not_ill_x(1:n_not_ill)=x0(ceil(2e3*rand(1,n_not_ill)));
M_not_ill_y(1:n_not_ill)=y0(ceil(2e3*rand(1,n_not_ill)));
%
%
% %%
% N_ill=0;
%
%
% %j=6;
% %old_people=Who_not_ill(Who_not_ill>groups(j-1) & Who_not_ill<=(groups(j)-N*0.15*0.1));
r_end=2e3;
% %ratio=Ratio_09(:,3);
for n=1:440
    if n>350
        n
    end
    
    r_end=ratio(n)*2e3;
    if n>200
        
        r_end=ratio(n)*1.15*2e3;
    end
    if n>290&&n<=350
        r_end=ratio(n)*1.2*2e3;
    end
    
    if n>350
        r_end=ratio(350)*2e3*1.1;
    end
    if n>380
        r_end=ratio(350)*2e3*1;
    end
    if n>395
        r_end=ratio(350)*2e3*1;
    end
    
    
    
    
    
    clear temp
    clear temp2
    temp=ceil(normrnd(14,4,1,length(Who_ill(:,1)))/2+...
        normrnd(18,6,1,length(Who_ill(:,1)))/2);%time to recover
    %a1=
    a1=Who_ill((n-Who_ill(:,2)-temp')>0,1);
    a2=Who_ill((n-Who_ill(:,2)-temp')>0,2);
    
    Who_rec(end+1:end+length(a1),1)=a1;
    Who_rec(end+1-length(a2):end,2)=a2;
    active_case=setdiff(Who_ill(:,1),Who_rec(:,1));%numer of active cases
    N_ill(n,1)=length(setdiff(Who_ill(:,1),Who_rec(:,1)));
    N_ill(n,2)=length(find(active_case>N*0.85));
    
    
    %add new sick people - chance of 20%
    if ratio(n)<1.3&&n<100
        ch=0.8;
    else
        ch=0.95;
    end
    if rand(1,1)>ch
        n;
        M_ill_x(length(M_ill_x)+1)=x0(ceil(r_end*rand(1,1)));
        M_ill_y(length(M_ill_y)+1)=y0(ceil(r_end*rand(1,1)));
        
        M_not_ill_x=M_not_ill_x(1:end-1);
        M_not_ill_y=M_not_ill_y(1:end-1);
        
        Who_ill(end+1,1)=Who_not_ill(ceil(rand(1,1)*length(Who_not_ill)));
        Who_ill(end,2)=n;
        Who_ill(end,3)=0;
        Who_not_ill=Who_not_ill(Who_not_ill~=Who_ill(end,1));
    end
    a =[255.1, 371.1, 4.477,36.17];
    if rand(1,1)<a(1)*exp(-(n-a(2)).^2/a(3)^2/12)/854
        1
        M_ill_x(length(M_ill_x)+1)=x0(ceil(r_end*rand(1,1)));
        M_ill_y(length(M_ill_y)+1)=y0(ceil(r_end*rand(1,1)));
        
        M_not_ill_x=M_not_ill_x(1:end-1);
        M_not_ill_y=M_not_ill_y(1:end-1);
        
        Who_ill(end+1,1)=Who_not_ill(ceil(rand(1,1)*length(Who_not_ill)));
        Who_ill(end,2)=n;
        Who_ill(end,3)=1;
        Who_not_ill=Who_not_ill(Who_not_ill~=Who_ill(end,1));
    end
     if n>380&&n<410
        a2 =[255.1, 410, 4.477,36.17];
        xx=a2(1)*exp(-(n-a2(2)).^2/a2(3)^2/12)/854/3;
        t1=find(Who_ill(:,2)==n-1);
        if xx>0
            n
            length(t1)
        end
        Who_ill(t1(rand(size(t1))<xx),3)=3;
    end
    
    %M_free_x=M_free_x(1:end-n_vaccine);
    %M_free_y=M_free_y(1:end-n_vaccine);
    
    %%
    l=3;
    n_not_ill=size(Who_not_ill,2);
    x1_not_ill=normrnd(0,2*0.5*(10^l),1,n_not_ill);
    y1_not_ill=normrnd(0,2*0.5*(10^l),1,n_not_ill);
    
    n_ill=size(Who_ill,1);
    
    x1_ill=normrnd(0,2*0.5*(10^l),1,n_ill);
    y1_ill=normrnd(0,2*0.5*(10^l),1,n_ill);
    
    temp=normrnd(5,1,1,length(Who_ill(:,2)));
    temp(temp<0)=0;
    
    temp2=n-Who_ill(:,2);%number of days from infection
    %temp3=find(rand(1,length(days))>0.25);%people with symptoms
    temp4=find(temp<temp2);%sick people that are still sick.
    %
    %
    %     %                 x1_ill(intersect(temp3,temp4))=0;
    %     %                 y1_ill(intersect(temp3,temp4))=0;
    %     %
    M_ill_x=M_ill_x+x1_ill;
    M_ill_y=M_ill_y+y1_ill;
    M_ill_x=mod(M_ill_x,r_end);M_ill_y=mod(M_ill_y,r_end);
    
    
    %%
    
    
    M_not_ill_x=M_not_ill_x+x1_not_ill;
    M_not_ill_y=M_not_ill_y+y1_not_ill;
    
    
    M_not_ill_x=mod(M_not_ill_x,r_end);
    M_not_ill_y=mod(M_not_ill_y,r_end);
    %print(M_free_y(M_free_y>2e3))
    %print(M_free_y(M_free_x>2e3))
    
    
    t_M_ill=length(M_ill_x);
    if t_M_ill>0
        M_ill_x1=M_ill_x'*ones(1,n_not_ill);
        M_ill_y1=M_ill_y'*ones(1,n_not_ill);
        
        
        if n_not_ill>0
            M_not_ill_x1=ones(n_ill,1)*M_not_ill_x;
            M_not_ill_y1=ones(n_ill,1)*M_not_ill_y;
            
            T_x=abs(M_not_ill_x1-M_ill_x1);
            T_y=abs(M_not_ill_y1-M_ill_y1);
            T=sqrt(T_x.^2+T_y.^2);
            %temp5=intersect(temp3,temp4);%people with symptoms after ~5 days
            
            map=ones(1,length(Who_ill));
            
            
            if ~isempty(Who_vac>0)
                
                map=ones(1,length(Who_ill));
                temp1=find(n-Who_vac(:,1)>10);
                [val,pos]=intersect(Who_ill,Who_vac(temp1,1));
                map(pos)=0.5;
            end
            
            
            x_n=n-Who_ill(:,2);
            %size(Who_ill)
            
            
            r_t_0=exp(-(x_n-7).^2/2)+exp(-(x_n-4).^2/2)+exp(-(x_n-5).^2/2)+exp(-(x_n-6).^2/2);%chance of being infection a function of the day.
            r_t_1=exp(-(x_n-2).^2/2)+exp(-(x_n-4).^2/2)+exp(-(x_n-5).^2/2)+exp(-(x_n-3).^2/2);%chance of being infection a function of the day.
            r_t=zeros(size(x_n));
            
            if size(Who_ill,2)>2
                r_t(Who_ill(:,3)==0)=r_t_0(Who_ill(:,3)==0);
                r_t(Who_ill(:,3)>0)=r_t_1(Who_ill(:,3)>0);
            else
                r_t=r_t_0;
            end
            r_t(r_t>1)=1;
            
            
            r_t(r_t>1)=1;
            
            risk=ones(size(r_t));
            if size(Who_ill,2)>2
                risk=Who_ill(:,3)*1+1;
            end
            %risk=ones(size(r_t));
            clear map
            map=ones(size(r_t));
            size(map);
            if ~isempty(Who_vac>0)
                
                temp1=find(n-(20+Who_vac(:,1))>10);
                [val,pos]=intersect(Who_ill(:,1),Who_vac(temp1,1));
                size(pos);
                map(pos)=0.5;
                %if n>380
                %   map(pos)=0.15;%Prevention of 40% infection of vaccinated
                %end
                
            end
            n;
            if size(r_t,1)~=size(gamma,1)
                gamma=gamma';
            end
            size(r_t);
            size(gamma);
            
            numeber_of_oc=ceil(rand(1,t_M_ill)*2);
            a_symptom=round(rand(1,t_M_ill))/2+0.5;%half of the people are asymptomatic with P=1/2;
            %                         n;
            T(T<3)=3; %keepind SD
            temp3=find(a_symptom==1);
            temp5=intersect(temp3,temp4);%people with symptoms after ~5 days
            
            %                         for jj=1:length(temp5)
            %                             a=T(temp5(jj),:);
            %                             a(a<8)=8;
            %                             T(temp5(jj),:)=a;
            %                             min(a);
            %                         end
            t=1;
            map=map(1:length(r_t));
            size(T);
            P1=exp(-T.^2/2/2.4^2).*r_t.*map*number_of_people.*numeber_of_oc'.*risk;
            P1_a=P1;
            if ~isempty(find((P1)>1e-3))
                [ii,jj]=find((P1)>1e-3);
                ii;
                P1_a=P1_a(ii,:);
                
                
                [M,I] = max(P1_a,[],1);
                size(I);
                J=risk(ii);
                J=J(I);
            else
                J=0;
            end
            
            %size(T)
            
            
            left_time1=n-(20+Who_vac(Who_not_ill,1));
            left_time2=n-(20+Who_vac(Who_not_ill,2));
            
            x=1-0.3.*(J>1);
            map1=(fittedmodel(left_time1)+.4*fittedmodel1(left_time2));
            if size(x,1)~=size(map1,1)
                1;
                x=x';
            end
            map=1-x.*map1;
            
            
            
            
            
            
            P=sum(0.7*P1_a,1);
            size(P);
            size(map);
            
            P=(P.*map');
            %save('P.mat','P','n','M_not_ill_x','map','x','map1')
            P=floor(P+rand(size(P)));
            temp=1*(J(P>0)>1);
            
            %P1=exp(-T.^2/2/2.4^2).*r_t.*map*number_of_people.*numeber_of_oc'.*risk;
            
            P_a=(sum(P1(risk==1,:),1)).*map'*0.7;
            P_b=(sum(P1(risk==2,:),1)).*map'*0.7;
            P_c=(sum(P1(risk>2,:),1)).*map'*0.7;
            
            P_a=floor(P_a+rand(size(P_a)));
            P_b=floor(P_b+rand(size(P_b)));
            P_c=floor(P_c+rand(size(P_c)));
            
            P_a(P_a>1)=1;
            P_b(P_b>1)=1;
            P_c(P_c>1)=1;
            
            P=P_a+P_b+P_c;
            P(P>1)=1;
            n;
           if n==395
               n
                save('P.mat','P_a','P_b','P_c','P','P1','risk','n')
           end
          
            
            
            new_M_not_ill_x=M_not_ill_x(P<1);
            new_M_not_ill_y=M_not_ill_y(P<1);
            j;
            sum(P);
           
            if sum(P>0)
                %                 if sum(temp>0)
                %                     save('J.mat','temp')
                %                 end
                M_ill_x(end+1:end+sum(P(P>0)))=M_not_ill_x(P>0);
                M_ill_y(end+1:end+sum(P(P>0)))=M_not_ill_y(P>0);
                Who_ill(end+1:end+sum(P(P>0)),1)=Who_not_ill(P>0);
                Who_ill(end+1-sum(P(P>0)):end,2)=n;
                p2=(P_a>0)*1+(P_b>0)*2+(P_c>0)*4;
                p2(p2>4)=4;
                p2=p2(p2>0)-1;
                Who_ill(end+1-sum(P(P>0)):end,3)=p2;
                Who_not_ill=setdiff(Who_not_ill,Who_ill(:,1));
            end
            n_not_ill=length(new_M_not_ill_x);
            M_not_ill_x=new_M_not_ill_x;
            M_not_ill_y=new_M_not_ill_y;
            %adding reinfection after 150 days:
            
            
        end
    end
    
    
end
% %end
%
%end