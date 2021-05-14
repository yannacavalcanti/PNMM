function [err_norm] = cmpt_err_norm(X,X0)

err_norm = norm(X-X0,'fro')/norm(X0,'fro');