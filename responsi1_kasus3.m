% =========================================================================
% RESPONSI 1 PRAKTIKUM METODE NUMERIK
% KASUS 3 : Persamaan Nonlinear - Bola di Dalam Tabung
% GNU Octave
%
%   Tabung menyinggung bola:  t = 2r
%   V bola   = (4/3) pi r^3
%   V tabung = pi r^2 t = 2 pi r^3
%   V bola : V tabung = (4/3) : 2 = 2 : 3
%
%   V bola = 300  ->  V tabung = 300 * 3/2 = 450
%   V tabung di luar bola = 450 - 300 = 150
%
%   V luar = (2/3) pi r^3 = 150
%   f(r)   = (2/3) pi r^3 - 150 = 0
% =========================================================================
1;   % penanda bahwa file ini adalah script (bukan file function)

% ------------------------- FUNGSI BANTU ----------------------------------
function [c, iter, tabel] = biseksi(f, a, b, tol, maxiter)
  % metode biseksi: c = (a+b)/2
  tabel = [];
  c = a;
  for iter = 1:maxiter
    clama = c;
    c  = (a + b) / 2;
    fa = f(a);  fb = f(b);  fc = f(c);
    if iter == 1
      ea = NaN;
    else
      ea = abs(c - clama);
    end
    tabel = [tabel; iter, a, b, c, fa, fb, fc, ea];
    if fa * fc < 0
      b = c;
    elseif fa * fc > 0
      a = c;
    else
      break;          % c tepat akar
    end
    if iter > 1 && ea < tol
      break;
    end
  end
end

function [c, iter, tabel] = regula_falsi(f, a, b, tol, maxiter)
  % metode regula-falsi: c = b - f(b)(b-a)/(f(b)-f(a))
  tabel = [];
  c = a;
  for iter = 1:maxiter
    clama = c;
    fa = f(a);  fb = f(b);
    c  = b - fb * (b - a) / (fb - fa);
    fc = f(c);
    if iter == 1
      ea = NaN;
    else
      ea = abs(c - clama);
    end
    tabel = [tabel; iter, a, b, c, fa, fb, fc, ea];
    if fa * fc < 0
      b = c;
    elseif fa * fc > 0
      a = c;
    else
      break;
    end
    if iter > 1 && ea < tol
      break;
    end
  end
end

function cetak_tabel(tabel)
  fprintf('   %4s %10s %10s %12s %12s %12s %12s %10s\n', ...
          'iter', 'a', 'b', 'c', 'f(a)', 'f(b)', 'f(c)', 'galat');
  for r = 1:size(tabel,1)
    if isnan(tabel(r,8))
      fprintf('   %4d %10.6f %10.6f %12.8f %12.5f %12.5f %12.6f %10s\n', ...
              tabel(r,1:7), '-');
    else
      fprintf('   %4d %10.6f %10.6f %12.8f %12.5f %12.5f %12.6f %10.2e\n', ...
              tabel(r,:));
    end
  end
end

% ------------------------- PROGRAM UTAMA ---------------------------------
clc; clear; close all;

% ---- 1) perumusan persamaan nonlinear ----
V_bola   = 300;                    % cm^3 (diketahui)
V_tabung = V_bola * 3/2;           % dari perbandingan 2 : 3
V_luar   = V_tabung - V_bola;      % volume tabung di luar bola

f = @(r) (2/3) * pi * r.^3 - V_luar;     % f(r) = (2/3) pi r^3 - 150

fprintf('==================================================================\n');
fprintf(' KASUS 3 : PERSAMAAN NONLINEAR (BOLA DI DALAM TABUNG)\n');
fprintf('==================================================================\n');
fprintf('\n 1) Perumusan persamaan\n');
fprintf('   V bola : V tabung = 4/3 pi r^3 : 2 pi r^3 = 2 : 3\n');
fprintf('   V bola              = %g cm^3\n', V_bola);
fprintf('   V tabung            = %g cm^3\n', V_tabung);
fprintf('   V tabung di luar bola = %g cm^3\n', V_luar);
fprintf('   Persamaan nonlinear : f(r) = (2/3) pi r^3 - %g = 0\n', V_luar);

% ---- 2) metode grafik untuk menentukan nilai awal ----
r = linspace(0, 8, 400);
figure(1);
plot(r, f(r), 'b-', 'LineWidth', 2); hold on;
plot([0 8], [0 0], 'k-');                          % sumbu f(r) = 0
a = 4;  b = 5;                                     % interval dari grafik
plot([a b], [f(a) f(b)], 'ro', 'MarkerFaceColor', 'r');
hold off; grid on;
xlabel('r (cm)'); ylabel('f(r)');
title('Metode grafik: f(r) = (2/3)\pi r^3 - 150');
legend('f(r)', 'f(r) = 0', 'a = 4 dan b = 5', 'Location', 'northwest');
print(gcf, 'kasus3_grafik.png', '-dpng', '-r120');

fprintf('\n 2) Metode grafik: grafik f(r) memotong sumbu r antara 4 dan 5\n');
fprintf('   a = %g -> f(a) = %.5f\n', a, f(a));
fprintf('   b = %g -> f(b) = %.5f\n', b, f(b));
fprintf('   f(a) * f(b) = %.5f < 0 -> akar berada di dalam [a, b]\n', f(a)*f(b));

% ---- 3) biseksi dan regula-falsi ----
tol     = 1e-6;
maxiter = 100;

fprintf('\n 3a) Metode Biseksi (tol = %g)\n\n', tol);
[r_bis, it_bis, tab_bis] = biseksi(f, a, b, tol, maxiter);
cetak_tabel(tab_bis);
fprintf('\n   Akar biseksi = %.8f cm  (%d iterasi)\n', r_bis, it_bis);

fprintf('\n 3b) Metode Regula-Falsi (tol = %g)\n\n', tol);
[r_rf, it_rf, tab_rf] = regula_falsi(f, a, b, tol, maxiter);
cetak_tabel(tab_rf);
fprintf('\n   Akar regula-falsi = %.8f cm  (%d iterasi)\n', r_rf, it_rf);

% ---- 4) hasil akhir dan verifikasi ----
r_akar = r_bis;
t_akar = 2 * r_akar;
r_eksak = (225/pi)^(1/3);          % penyelesaian analitik r^3 = 225/pi
fprintf('\n 4) Hasil akhir (memakai hasil biseksi)\n');
fprintf('   jari-jari bola & tabung r = %.6f cm\n', r_akar);
fprintf('   tinggi tabung t = 2r      = %.6f cm\n', t_akar);
fprintf('   V bola                    = %.6f cm^3\n', 4/3*pi*r_akar^3);
fprintf('   V tabung                  = %.6f cm^3\n', pi*r_akar^2*t_akar);
fprintf('   V tabung di luar bola     = %.6f cm^3\n', pi*r_akar^2*t_akar - 4/3*pi*r_akar^3);
fprintf('   pembanding analitik r = (225/pi)^(1/3) = %.8f cm\n', r_eksak);
fprintf('==================================================================\n');
