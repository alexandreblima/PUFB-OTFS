function [mse,frob,sparsity,full,stability]=r22d_test(methods,p,physical,output)
    directory=fullfile(output,'test_channels');if ~isfolder(directory),mkdir(directory);end
    manifest_path=fullfile(output,'r22d_test_manifest.csv');manifest=r22d_model('read',manifest_path);
    assert(height(manifest)==1000&&isequal(manifest.Seed,p.seeds.TEST+1000*repelem((1:5)',200)+repmat((1:200)',5,1)));
    assert(isequal(manifest.VelocityIndex,repelem((1:5)',200))&&isequal(manifest.ChannelID,repmat((1:200)',5,1)));
    mse_rows=cell(48000,8);full_rows=cell(48000,5);frob_rows=cell(6000,5);sparse_rows=cell(6000,11);
    metric_cursor=0;channel_cursor=0;
    for index=1:height(manifest)
        path=fullfile(directory,sprintf('channel_%04d.mat',index));
        if isfile(path)
            loaded=load(path,'record','seed');record=loaded.record;assert(loaded.seed==manifest.Seed(index));
        else
            seed=manifest.Seed(index);velocity=manifest.VelocityIndex(index);G=generate_channel_r20a(24,12,p.velocity_mps(velocity),seed,physical);
            record=struct('SparseMSE',zeros(6,8),'FullMSE',zeros(6,8),'Frobenius',zeros(6,1), ...
                'Diagnostics',zeros(6,7),'OracleMasks',false(288^2,6));
            for method=1:6
                U=methods(method).U;H=U'*G*U;mask=methods(method).Mask;
                assert(nnz(mask)==p.K,'R22D:TestK','TEST support changed.');
                [record.SparseMSE(method,:),record.FullMSE(method,:)]=confirm_exact_mse(H,mask,p.snr_db);
                record.Frobenius(method)=norm(H-H.*mask,'fro')/norm(H,'fro');
                energy=abs(H(:)).^2;probability=energy/sum(energy);[ordered,order]=sort(probability,'descend');cumulative=cumsum(ordered);
                oracle=false(288^2,1);oracle(order(1:p.K))=true;record.OracleMasks(:,method)=oracle;
                coverage=nnz(oracle&mask(:))/p.K;
                entropy=-sum(probability(probability>0).*log(probability(probability>0)))/log(numel(probability));
                magnitude=abs(H(:));hoyer=(sqrt(numel(magnitude))-sum(magnitude)/norm(magnitude))/(sqrt(numel(magnitude))-1);
                record.Diagnostics(method,:)=[find(cumulative>=.95,1),find(cumulative>=.99,1),find(cumulative>=.999,1), ...
                    sum(probability.^2),entropy,hoyer,coverage];
            end
            assert(max(abs(record.FullMSE-record.FullMSE(1,:)),[],'all')<p.full_invariance_tolerance);
            save(path,'record','seed','G','-v7.3');
        end
        velocity=manifest.VelocityIndex(index);channel=manifest.ChannelID(index);manifest.ChannelGenerated(index)=true;
        for method=1:6
            channel_cursor=channel_cursor+1;name=methods(method).Name;
            frob_rows(channel_cursor,:)={name,velocity,channel,p.K,record.Frobenius(method)};
            sparse_rows(channel_cursor,:)=[{name,velocity,channel,p.K},num2cell(record.Diagnostics(method,:))];
            for point=1:8
                metric_cursor=metric_cursor+1;
                mse_rows(metric_cursor,:)={name,velocity,channel,p.snr_db(point),p.K,record.SparseMSE(method,point), ...
                    record.FullMSE(method,point),record.SparseMSE(method,point)-record.FullMSE(method,point)};
                full_rows(metric_cursor,:)={name,velocity,channel,p.snr_db(point),record.FullMSE(method,point)};
            end
        end
        if mod(index,10)==0||index==height(manifest)
            writetable(manifest,manifest_path);fprintf('R22D TEST exact solves %d/1000\n',index);
        end
    end
    mse=cell2table(mse_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','RetainedEntries','MSESparse','MSEFull','DeltaMSE'});
    full=cell2table(full_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','SNRdB','FullMSE'});
    frob=cell2table(frob_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','FrobeniusError'});
    sparsity=cell2table(sparse_rows,'VariableNames',{'Method','VelocityIndex','ChannelID','RetainedEntries','K95','K99','K999', ...
        'J4','NormalizedEntropy','Hoyer','SupportCoverage'});
    writetable(mse,fullfile(output,'r22d_test_exact_mse.csv'));writetable(full,fullfile(output,'r22d_full_receiver_mse.csv'));
    writetable(frob,fullfile(output,'r22d_test_frobenius.csv'));writetable(sparsity,fullfile(output,'r22d_test_sparsity.csv'));
    stability=support_stability(methods,p,output);writetable(stability,fullfile(output,'r22d_test_support_stability.csv'));
end

function T=support_stability(methods,p,output)
    count=40;rows=cell(30,6);cursor=0;
    for method=1:6
        masks=cell(5,1);
        for velocity=1:5
            masks{velocity}=false(288^2,count);
            for channel=1:count
                index=(velocity-1)*p.test_per_velocity+channel;
                loaded=load(fullfile(output,'test_channels',sprintf('channel_%04d.mat',index)),'record');
                masks{velocity}(:,channel)=loaded.record.OracleMasks(:,method);
            end
        end
        for velocity=1:5
            current=double(masks{velocity});intersection=current'*current;jaccard=intersection./(2*p.K-intersection);
            within=jaccard(triu(true(count),1));cross=[];
            for other=1:5
                if other==velocity,continue;end
                intersection=current'*double(masks{other});values=intersection./(2*p.K-intersection);
                cross=[cross;values(:)]; %#ok<AGROW>
            end
            cursor=cursor+1;rows(cursor,:)={methods(method).Name,velocity,mean(within),mean(cross),count,p.K};
        end
    end
    T=cell2table(rows,'VariableNames',{'Method','VelocityIndex','WithinVelocityJaccard','CrossVelocityJaccard','DiagnosticChannelsPerVelocity','RetainedEntries'});
end