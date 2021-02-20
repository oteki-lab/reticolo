% ===== Info =====
% Function for Send Mail
% URL: https://www.umi-mori.jp/article/matlab/send_mail.php
% Github: https://github.com/Masumi-M/SendMail_matlab

%body = "reticolo Simulation Done.";
%sendMail(body);

% ===== Main =====
function [] = sendMail(body,attachmentFilePath)
    % Modify these two lines to reflect your account and password.
    email = '';
    password = '';
    recipients = '';

    setpref('Internet', 'E_mail', email);
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'SMTP_Username', email);
    setpref('Internet', 'SMTP_Password', password);

    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');  %Socket for SSL (465)

    subject = 'Matlab Notification.';
    if attachmentFilePath == ""
        sendmail(recipients, subject, body);
    else
        sendmail(recipients, subject, body, attachmentFilePath);
    end

    return

end