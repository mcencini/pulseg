function tests = test_v22_fields
% TEST_V22_FIELDS The optional instance fields of PulSeg 2.2.0 in validate_ir and write_seq.
%
% Requires the MATLAB Pulseq toolbox on the path.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    testCase.assumeTrue( ...
        exist('mr.Sequence', 'class') == 8 || exist('mr.Sequence', 'file') == 2, ...
        'Pulseq MATLAB toolbox was not found on the MATLAB path.');
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.TestData.tmpDir = tmpDir;
    seqFile = fullfile(tmpDir, 'demo.seq');
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'demo'));
    create_demo_sequence(seqFile);
    testCase.TestData.ir = pulseg.import(seqFile, 'verbose', false);
end


function teardownOnce(testCase)
    if isfolder(testCase.TestData.tmpDir)
        rmdir(testCase.TestData.tmpDir, 's');
    end
end


function testImportWritesVersion220(testCase)
    testCase.verifyEqual(testCase.TestData.ir.pulseg_version, '2.2.0');
    pulseg.validate_ir(testCase.TestData.ir);
end


function testFlagsNeedVersion220(testCase)
    ir = testCase.TestData.ir;
    ir.pulseg_version = '2.1-alpha';
    ir.execution_stream(1).norot = 1;
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testFlagMustBeZeroOrOne(testCase)
    ir = testCase.TestData.ir;
    ir.execution_stream(1).nopos = 2;
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testNorotIsWrittenAsALabelWhereTheInstanceSetsIt(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        ir.execution_stream(k).norot = double(k == numel(ir.execution_stream));
    end
    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'norot.seq'));
    labels = seq.evalLabels('evolution', 'blocks');
    testCase.verifyTrue(isfield(labels, 'NOROT'));
    testCase.verifyEqual(labels.NOROT(1), 0);
    testCase.verifyEqual(labels.NOROT(end), 1);
end


function testRfShimIsWrittenWithEachRfBlock(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        n_rf = numel(ir.execution_stream(k).rf_amplitude);
        ir.execution_stream(k).rf_shim = repmat([1, 0.5i], n_rf, 1);
    end
    pulseg.validate_ir(ir);
    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'shim.seq'));
    testCase.verifyGreaterThan(seq.duration(), 0);
end


function testReadoutPlayedButNotStoredIsRefused(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        n_adc = numel(ir.execution_stream(k).adc_phase_offset);
        if n_adc > 0
            ir.execution_stream(k).adc_acquire = zeros(n_adc, 1);
            break;
        end
    end
    testCase.verifyError(@() pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'adc.seq')), ?MException);
end


function testShotIndexBeyondTheVariantsIsRefused(testCase)
    ir = with_two_shot_variants(testCase.TestData.ir);
    [k, idx] = first_use_of(ir, ir.base_blocks(arbitrary_gx_block(ir)).id);
    ir.execution_stream(k).gradient_shot_index(idx) = 2;
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testUnnormalizedShotVariantIsRefused(testCase)
    ir = with_two_shot_variants(testCase.TestData.ir);
    p = arbitrary_gx_block(ir);
    ir.base_blocks(p).block.gx.shot_waveforms(2, :) = 2 * ir.base_blocks(p).block.gx.shot_waveforms(2, :);
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testEachGradientPlaysTheShotItsInstanceSelects(testCase)
    ir = with_two_shot_variants(testCase.TestData.ir);
    p = arbitrary_gx_block(ir);
    [k, idx, j] = first_use_of(ir, ir.base_blocks(p).id);
    ir.execution_stream(k).gradient_shot_index(idx) = 1;
    pulseg.validate_ir(ir);

    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'shots.seq'));
    gx = seq.getBlock(seq_block_number(ir, k, j)).gx;
    variant = ir.base_blocks(p).block.gx.shot_waveforms(2, :);
    expected = ir.execution_stream(k).gradient_amplitude(idx, 1) * variant;
    testCase.verifyEqual(gx.waveform(:).', expected, 'AbsTol', 1e-4 * max(abs(expected)));
end


function testRotationMatrixIsWrittenAsARotationExtension(testCase)
    ir = testCase.TestData.ir;
    p = arbitrary_gx_block(ir);
    [k, idx, j] = first_use_of(ir, ir.base_blocks(p).id);
    n_grad = size(ir.execution_stream(k).gradient_amplitude, 1);
    R = [0 -1 0; 1 0 0; 0 0 1];
    ir.execution_stream(k).rotation_matrix = repmat(eye(3), 1, 1, n_grad);
    ir.execution_stream(k).rotation_matrix(:, :, idx) = R;
    pulseg.validate_ir(ir);

    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'rotation.seq'));
    block = seq.getBlock(seq_block_number(ir, k, j));
    testCase.verifyTrue(isfield(block, 'rotation'));
    expected = mr.makeRotation(R).rotQuaternion;
    testCase.verifyEqual(abs(dot(block.rotation.rotQuaternion, expected)), 1, 'AbsTol', 1e-9);
end


function p = arbitrary_gx_block(ir)
    for p = 1:numel(ir.base_blocks)
        b = ir.base_blocks(p).block;
        if isfield(b, 'gx') && ~isempty(b.gx) && strcmp(b.gx.type, 'grad') ...
                && isfield(b.gx, 'waveform') && numel(b.gx.waveform) > 2
            return;
        end
    end
    error('The demo IR has no base block with an arbitrary gx gradient.');
end


function ir = with_two_shot_variants(ir)
    p = arbitrary_gx_block(ir);
    g = ir.base_blocks(p).block.gx;
    w = g.waveform(:).';
    ramp = linspace(0.5, 1, numel(w));
    peak = max(abs(w .* ramp));
    g.shot_waveforms = [w; w .* ramp / peak];
    g.shot_first = [g.first, 0.5 * g.first / peak];
    g.shot_last = [g.last, g.last / peak];
    ir.base_blocks(p).block.gx = g;
    for k = 1:numel(ir.execution_stream)
        n_grad = size(ir.execution_stream(k).gradient_amplitude, 1);
        ir.execution_stream(k).gradient_shot_index = zeros(n_grad, 1);
    end
end


function [k, idx, j] = first_use_of(ir, base_id)
% FIRST_USE_OF First instance using a base block, with its gradient-event index and block position.
    vs_ids = [ir.virtual_segments.id];
    base_ids = [ir.base_blocks.id];
    for k = 1:numel(ir.execution_stream)
        inst = ir.execution_stream(k);
        vs = ir.virtual_segments(vs_ids == inst.virtual_segment_id);
        idx = 0;
        for j = 1:numel(vs.base_block_ids)
            id = vs.base_block_ids(j);
            if id < 2
                continue;
            end
            b = ir.base_blocks(base_ids == id).block;
            if any(cellfun(@(a) isfield(b, a) && ~isempty(b.(a)), {'gx', 'gy', 'gz'}))
                idx = idx + 1;
                if id == base_id && inst.gradient_amplitude(idx, 1) ~= 0
                    return;
                end
            end
        end
    end
    error('No instance plays base block %d with a non-zero gx amplitude.', base_id);
end


function n = seq_block_number(ir, k, j)
    vs_ids = [ir.virtual_segments.id];
    n = j;
    for m = 1:k-1
        vs = ir.virtual_segments(vs_ids == ir.execution_stream(m).virtual_segment_id);
        n = n + numel(vs.base_block_ids);
    end
end
