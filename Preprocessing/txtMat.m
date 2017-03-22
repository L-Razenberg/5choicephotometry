function txtMat()

lst=dir('*T.txt')
for i=1:numel(lst)
    fname=lst(i).name;
    if fname(end-4)=='T'
        data=load(fname);
        data=data(:,1:3);
        eval(['save ' fname(1:end-3) 'mat data'])
    end
end

end