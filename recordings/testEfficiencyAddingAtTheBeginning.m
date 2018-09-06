%https://stackoverflow.com/questions/40035629/add-a-new-element-to-the-beginning-of-an-existing-cell-array
%Compare which way of adding an element at the beginning of an array is
%faster:
a = [0] ; b = [0]; c = [0];
t1 = [0]; t2 = [0]; t3 = [0];
for i = 1:10000
    disp(['i = ',num2str(i)])
    tic
        a = [i a];
    t1(end+1) = toc;
    
    tic
%         b = fliplr(b) ; b(end+1) = i ; b = fliplr(b);
        b = testAddingAtTheBeginning(b,i);
    t2(end+1) = toc;
    
    tic
        c = horzcat(i,c);
    t3(end+1) = toc;
end
figure ; plot(a,t1,'-',a,t2,'-',a,t3,'-')
legend('vec = [elementToAdd vec]','flipvecMyFunction then vec(end+1) then flipvecMyFunction', ...
       'vec = horzcat(elementToAdd,vec)')