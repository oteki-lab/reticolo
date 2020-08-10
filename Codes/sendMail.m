% ===== Info =====
% Function for Send Mail
% URL: https://www.umi-mori.jp/article/matlab/send_mail.php
% Github: https://github.com/Masumi-M/SendMail_matlab

%body = "reticolo Simulation Done.";
%sendMail(body);

% ===== Main =====
function [] = sendMail(body)
    % Modify these two lines to reflect your account and password.
    email = '';
    password = '';
    recipients = email;

    setpref('Internet', 'E_mail', email);
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'SMTP_Username', email);
    setpref('Internet', 'SMTP_Password', password);

    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');  %Socket for SSL (465)

    subject = 'Matlab Notification.';
    message_line1 = strcat('Matlab PGM Ended.');
    message_line2 = strcat('Attachment: [ Att ]');
    message_line3 = strcat('Datetime: [ ', datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), ' ]');

    sendmail(recipients, subject, body);

    return

end