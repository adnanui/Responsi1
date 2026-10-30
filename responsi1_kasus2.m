% =========================================================================
% RESPONSI 1 PRAKTIKUM METODE NUMERIK
% KASUS 2 : SPL - Harga Bahan Makanan per Kilogram
% GNU Octave
%
%   Variabel : x1 = harga beras/kg, x2 = harga telur/kg,
%              x3 = harga tepung/kg, x4 = harga ayam/kg
%   Metode   : cek dominan diagonal -> Gauss-Seidel (metode tidak langsung)
% =========================================================================
1;   % penanda bahwa file ini adalah script (bukan file function)

% ------------------------- FUNGSI BANTU ----------------------------------
function [x, iter, hist] = gauss_seidel(A, b, x0, tol, maxiter)
  % Gauss-Seidel: x_i = ( b_i - sum_{j~=i} a_ij x_j ) / a_ii
  % berhenti jika galat relatif semu maksimum < tol atau iter = maxiter
  n    = length(b);
  x    = x0;
  hist = [];
  for iter = 1:maxiter
    xlama = x;
    for i = 1:n
      j    = [1:i-1, i+1:n];
      x(i) = (b(i) - A(i,j) * x(j)) / A(i,i);
    end
    ea   = max(abs((x - xlama) ./ x));
    hist = [hist; iter, x', ea];
    if ea < tol
      break;
    end
  end
end

% ------------------------- PROGRAM UTAMA ---------------------------------
clc; clear; close all;

% ---- 1) data soal, semua satuan dikonversi ke kg ----
% konversi: 1 kg = 10 hg = 100 dag = 1000 g
%            beras     telur      tepung     ayam
data = [ 2.5,       1500/1000, 5,         500/1000;    % Adel   (kg, g, kg, g)
         40/10,     0.75,      1000/1000, 2;           % Briza  (hg, kg, g, kg)
         2000/1000, 15/10,     50/100,    6000/1000;   % Caca   (g, hg, dag, g)
         1,         4250/1000, 20/10,     9/10 ];      % Dinda  (kg, g, hg, hg)
total = [161000; 149000; 251500; 173000];
nama  = {'Adel', 'Briza', 'Caca', 'Dinda'};
bahan = {'Beras', 'Telur', 'Tepung', 'Ayam'};

A = data;
b = total;
n = length(b);

fprintf('==================================================================\n');
fprintf(' KASUS 2 : SPL HARGA BAHAN MAKANAN\n');
fprintf('==================================================================\n');
fprintf('\n 1) Bentuk SPL (satuan kg; x1=beras, x2=telur, x3=tepung, x4=ayam)\n\n');
for i = 1:n
  fprintf('   %-6s: %5.2f x1 + %5.2f x2 + %5.2f x3 + %5.2f x4 = %9.0f\n', ...
          nama{i}, A(i,1), A(i,2), A(i,3), A(i,4), b(i));
end

fprintf('\n   Bentuk matriks A x = b\n');
disp('   A ='); disp(A);
disp('   b ='); disp(b);

% ---- 2) analisis syarat kekonvergenan Gauss-Seidel ----
% syarat cukup: matriks dominan diagonal secara ketat
%               |a_ii| > sum_{j~=i} |a_ij|   untuk setiap baris i
fprintf('\n 2) Analisis dominan diagonal (urutan persamaan awal)\n');
for i = 1:n
  jml = sum(abs(A(i,:))) - abs(A(i,i));
  if abs(A(i,i)) > jml
    ket = 'dominan';
  else
    ket = 'TIDAK dominan';
  end
  fprintf('   baris %d (%s): |a%d%d| = %.2f  vs  jumlah lainnya = %.2f -> %s\n', ...
          i, nama{i}, i, i, abs(A(i,i)), jml, ket);
end

% cari penataan ulang baris agar tiap persamaan dominan di satu kolom berbeda
kol_dom = zeros(n,1);
for i = 1:n
  for j = 1:n
    if abs(A(i,j)) > sum(abs(A(i,:))) - abs(A(i,j))
      kol_dom(i) = j;
    end
  end
end

if all(kol_dom > 0) && numel(unique(kol_dom)) == n
  urut = zeros(1, n);
  urut(kol_dom) = 1:n;          % persamaan ke-urut(j) ditaruh di baris j
  Ab = A(urut, :);
  bb = b(urut);
  fprintf('\n   Matriks dapat dibuat dominan diagonal dengan menukar urutan persamaan:\n');
  for j = 1:n
    fprintf('   persamaan baris %d <- %s (dominan pada kolom %s)\n', ...
            j, nama{urut(j)}, bahan{j});
  end
  disp('   A (setelah ditukar) ='); disp(Ab);
  disp('   b (setelah ditukar) ='); disp(bb);
  bisa_iteratif = true;
else
  Ab = A;  bb = b;
  bisa_iteratif = false;
  fprintf('\n   Tidak ada penataan baris yang membuat matriks dominan diagonal.\n');
end

% syarat perlu dan cukup: jari-jari spektral matriks iterasi Gauss-Seidel < 1
function rho = rho_gs(M)
  D = diag(diag(M));  L = tril(M, -1);  U = triu(M, 1);
  T = -((D + L) \ U);
  rho = max(abs(eig(T)));
end
fprintf('\n   Jari-jari spektral matriks iterasi Gauss-Seidel:\n');
fprintf('     urutan awal         : %.4f\n', rho_gs(A));
fprintf('     setelah ditukar     : %.4f  (< 1 berarti konvergen)\n', rho_gs(Ab));

% ---- 3) penyelesaian ----
tol     = 1e-6;     % toleransi galat relatif semu
maxiter = 200;
x0      = zeros(n,1);

if bisa_iteratif
  fprintf('\n 3) Penyelesaian dengan Gauss-Seidel (x0 = 0, tol = %g, maks iter = %d)\n\n', ...
          tol, maxiter);
  [x, iter, hist] = gauss_seidel(Ab, bb, x0, tol, maxiter);

  fprintf('   %4s %14s %14s %14s %14s %12s\n', 'iter', 'x1 (beras)', ...
          'x2 (telur)', 'x3 (tepung)', 'x4 (ayam)', 'galat');
  for r = 1:size(hist,1)
    fprintf('   %4d %14.4f %14.4f %14.4f %14.4f %12.2e\n', hist(r,:));
  end
  fprintf('\n   Konvergen pada iterasi ke-%d\n', iter);
else
  fprintf('\n 3) Penyelesaian dengan eliminasi Gauss (metode langsung)\n');
  x = Ab \ bb;     % bisa diganti dengan program eliminasi Gauss / Gauss-Jordan
end

% ---- 4) hasil akhir dan verifikasi ----
x_langsung = A \ b;      % pembanding (metode langsung bawaan Octave)
fprintf('\n 4) Hasil akhir: harga bahan makanan per kilogram\n');
for j = 1:n
  fprintf('   %-7s : Rp %12.2f /kg   (pembanding A\\b: Rp %.2f)\n', ...
          bahan{j}, x(j), x_langsung(j));
end
fprintf('\n   Cek residu max|A x - b| = %.3e\n', max(abs(A*x - b)));
fprintf('==================================================================\n');
