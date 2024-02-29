
%�˼�����λ��
rod_point = [rod_x, 0.5.*rod_h];
%��Ƥ��Ԫ�ĳ���
for ii=1:length(rod_point)-1
    len_b(ii,1) = norm(rod_point(ii+1,:)-rod_point(ii,:));
end
%���嵥Ԫ�ĳ���
for ii=1:length(web_loc)
    len_w(ii,1) = rod_h(web_loc(ii));
end
%�����������
room_area = zeros(length(web_loc)-1, 1);
for ii=1:length(web_loc)-1
    for jj = web_loc(ii):web_loc(ii+1)-1
        room_area(ii) = room_area(ii)+ ...
            0.5 * (len_w(ii)+len_w(ii+1)) * (rod_x(jj+1)-rod_x(jj));
    end
end

%�����������Ծ�
Jx = 0;
for ii = 1:length(rod_h)
    Jx = Jx + 2 * rod_A(ii) * (0.5*rod_h(ii))^2;
end

%ͨ��ÿ���˼���ļ�������
delta_q_b = -(Q_all/Jx).* rod_A .* (0.5.*rod_h);
q_b0 = zeros(length(rod_x)-1, 1);
q_b0(1) = delta_q_b(1);
%���踹��û�м���
for ii=2:length(rod_x)
    q_b0(ii) = q_b0(ii-1) + delta_q_b(ii);
end

A_bend = zeros(length(web_loc)-1, length(web_loc)-1);
b_bend = zeros(length(web_loc)-1, 1);
%��ÿ�������е�λ�غɷ���
for ii=1:length(web_loc)-1
    b_num = web_loc(ii):web_loc(ii+1);%��Ƥ��Ԫ�����
    for jj=1:length(b_num)-1
        A_bend(ii,1:ii) = A_bend(ii,1:ii) + 2*len_b(b_num(jj))/t_b(b_num(jj));
    end
    A_bend(ii,ii) = A_bend(ii,ii) + len_w(ii)/t_w(ii);
    if ii~=length(web_loc)-1
        A_bend(ii,ii+1) = A_bend(ii,ii+1) - len_w(ii+1)/t_w(ii+1);
    end
    b_bend(ii) = 0;
    for jj=1:length(b_num)-1
        b_bend(ii) = b_bend(ii) + q_b0(b_num(jj)) * 2*len_b(b_num(jj))/t_b(b_num(jj));
    end
    
    %���һ�����ҵ�b���⴦��
    if ii == length(web_loc)-1
        b_bend(ii) = b_bend(ii) + q_b0(end) * len_w(ii+1)/t_w(ii+1);
        A_bend(ii,:) = A_bend(ii,:) + len_w(ii+1)/t_w(ii+1);
    end
end
%����ڼ��������µĸ������
q_w_bend = A_bend\b_bend;

%������Ƥ��Ԫ�ļ���
q_w_int = 0;
for ii=1:length(web_loc)-1
    q_w_int = q_w_int + q_w_bend(ii);
    for jj=web_loc(ii):web_loc(ii+1)-1
        q_b_bend(jj,1) = q_b0(jj) - q_w_int;
    end
end
q_w_bend = [q_w_bend; q_b0(end)-q_w_int];

%�����������
M_all = 0;
F_all = 0;
for ii=1:length(rod_x)-1
    M_all = M_all - 2 * q_b_bend(ii)...
        * norm(rod_point(ii+1,:)-rod_point(ii,:))...
        * distance(rod_point(ii,:), rod_point(ii+1,:), [0,0]);
    F_all = F_all + 2 * q_b_bend(ii) * (rod_point(ii+1,2)-rod_point(ii,2));          
end
for ii=1:length(web_loc)
    M_all = M_all - rod_x(web_loc(ii)) *  rod_h(web_loc(ii)) * q_w_bend(ii);
    F_all = F_all - rod_h(web_loc(ii)) * q_w_bend(ii);
end
xs = M_all/F_all

phi=1;%��Ťת��Ϊ1
A_torq = zeros(length(web_loc)-1, length(web_loc)-1);
b_torq = zeros(length(web_loc)-1, 1);
%��ÿ�������е�λ�غɷ���
for ii=1:length(web_loc)-1
    b_num = web_loc(ii):web_loc(ii+1);%��Ƥ��Ԫ�����
    for jj=1:length(b_num)-1
        A_torq(ii,1:ii) = A_torq(ii,1:ii) + 2*len_b(b_num(jj))/t_b(b_num(jj));
    end
    A_torq(ii,ii) = A_torq(ii,ii) + len_w(ii)/t_w(ii);
    if ii~=length(web_loc)-1
        A_torq(ii,ii+1) = A_torq(ii,ii+1) - len_w(ii+1)/t_w(ii+1);
    end
    b_torq(ii) = -phi*G*2*room_area(ii);
       
    %���һ�����ҵ����⴦��
    if ii == length(web_loc)-1
        A_torq(ii,:) = A_torq(ii,:) + len_w(ii+1)/t_w(ii+1);
    end
end
q_w_torq = A_torq\b_torq;

%Ť�ز�������Ƥ����
q_w_int = 0;
for ii=1:length(web_loc)-1
    q_w_int = q_w_int + q_w_torq(ii);
    for jj=web_loc(ii):web_loc(ii+1)-1
        q_b_torq(jj,1) =  - q_w_int;
    end
end
q_w_torq = [q_w_torq; -q_w_int];

%Ťת�նȼ���
M_all = 0;
F_all = 0;
for ii=1:length(rod_x)-1
    M_all = M_all + 2 * q_b_torq(ii)...
        * norm(rod_point(ii+1,:)-rod_point(ii,:))...
        * distance(rod_point(ii,:), rod_point(ii+1,:), [0,0]);
    F_all = F_all + 2 * q_b_torq(ii) * (rod_point(ii+1,2)-rod_point(ii,2));          
end
for ii=1:length(web_loc)
    M_all = M_all + rod_x(web_loc(ii)) *  rod_h(web_loc(ii)) * q_w_torq(ii);
    F_all = F_all - rod_h(web_loc(ii)) * q_w_torq(ii);
end
GJ = M_all/phi;

%Ť�ز�������Ƥ�븹���ʵ�ʼ���
q_b_torq = q_b_torq./M_all;
q_w_torq = q_w_torq./M_all;

%������飺
%ÿ�����ҵ�Ťת�ǣ������Ƿ����
% delta_phi(:,1)��Ť�������Ťת��
% delta_phi(:,2)�Ǽ��������Ťת��
% ��Ť�� Torque_check 
% �ܼ��� Q_check
for ii=1:length(web_loc)-1
    delta_phi(ii,1) = 0;
    delta_phi(ii,2) = 0;
    for jj=web_loc(ii):web_loc(ii+1)-1
        delta_phi(ii,1) = delta_phi(ii,1) + 2*q_b_torq(jj)*len_b(jj)/t_b(jj);
        delta_phi(ii,2) = delta_phi(ii,2) + 2*q_b_bend(jj)*len_b(jj)/t_b(jj);
    end
    delta_phi(ii,1) = delta_phi(ii,1) - q_w_torq(ii)*len_w(ii)/t_w(ii);
    delta_phi(ii,1) = delta_phi(ii,1) + q_w_torq(ii+1)*len_w(ii+1)/t_w(ii+1);
    delta_phi(ii,2) = delta_phi(ii,2) - q_w_bend(ii)*len_w(ii)/t_w(ii);
    delta_phi(ii,2) = delta_phi(ii,2) + q_w_bend(ii+1)*len_w(ii+1)/t_w(ii+1);
    delta_phi(ii,:) = delta_phi(ii,:)./(2*room_area(ii))./G;
end     
Q_check = 0;
Torque_check = 0;
q_b = q_b_bend + q_b_torq;
q_w = q_w_bend + q_w_torq;
for ii=1:length(rod_x)-1
    Torque_check = Torque_check + 2 * q_b(ii)...
        * norm(rod_point(ii+1,:)-rod_point(ii,:))...
        * distance(rod_point(ii,:), rod_point(ii+1,:), [xs,0]);
    Q_check = Q_check + 2 * q_b(ii) * (rod_point(ii+1,2)-rod_point(ii,2));          
end
for ii=1:length(web_loc)
    Torque_check = Torque_check + (rod_x(web_loc(ii))-xs) *  rod_h(web_loc(ii)) * q_w(ii);
    Q_check = Q_check - rod_h(web_loc(ii)) * q_w(ii);
end
% Torque_check/delta_phi(1,1)

figure(1);
title('wing box');hold on;
plot(rod_x, 0.5.*rod_h, '-o', 'LineWidth', 4);hold on;
plot(rod_x, -0.5.*rod_h, '-o', 'LineWidth', 4);hold on;
for ii=1:length(web_loc)
    plot([rod_x(web_loc(ii)), rod_x(web_loc(ii))], 0.5.*[rod_h(web_loc(ii)),...
        -rod_h(web_loc(ii))], '-', 'LineWidth', 4);
    hold on;
end
axis equal;

    

    
    
    
    
    
    
    
