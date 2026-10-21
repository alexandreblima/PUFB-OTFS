function result=r22a_io(action,path,value)
    if nargin<3,value=[];end
    result=[];
    switch action
        case 'json'
            fid=fopen(path,'w');assert(fid>=0,'R22A:Write','Cannot write output.');cleanup=onCleanup(@()fclose(fid));
            fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
        case 'text'
            fid=fopen(path,'w');assert(fid>=0,'R22A:Write','Cannot write output.');cleanup=onCleanup(@()fclose(fid));
            fprintf(fid,'%s\n',value);
        case 'hash'
            fid=fopen(path,'rb');assert(fid>=0,'R22A:Hash','Cannot read input.');cleanup=onCleanup(@()fclose(fid));
            bytes=fread(fid,Inf,'*uint8');digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
            result=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
        case 'verify'
            entries=readtable(path,'Delimiter',',','VariableNamingRule','preserve','TextType','string');
            for file_index=1:height(entries)
                assert(strcmp(r22a_io('hash',char(entries.AbsolutePath(file_index))),char(entries.SHA256(file_index))), ...
                    'R22A:InputChanged','Hash changed: %s',entries.AbsolutePath(file_index));
            end
        case 'manifest'
            files=[dir(fullfile(value,'*'));dir(fullfile(path,'**','*'))];rows=cell(0,3);
            for file_index=1:numel(files)
                file=files(file_index);
                if file.isdir||any(strcmp(file.name,{'artifact_manifest.csv','artifact_manifest.json'})),continue;end
                absolute=fullfile(file.folder,file.name);rows(end+1,:)={string(absolute),string(r22a_io('hash',absolute)),file.bytes}; %#ok<AGROW>
            end
            result=cell2table(rows,'VariableNames',{'AbsolutePath','SHA256','Bytes'});
            writetable(result,fullfile(path,'artifact_manifest.csv'));
            r22a_io('json',fullfile(path,'artifact_manifest.json'),struct('entries',height(result), ...
                'manifest_csv_sha256',r22a_io('hash',fullfile(path,'artifact_manifest.csv'))));
        otherwise
            error('R22A:IOAction','Unknown IO operation.');
    end
end