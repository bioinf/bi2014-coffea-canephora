



<!DOCTYPE html>
<html lang="en" class="">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta http-equiv="Content-Language" content="en">
    
    
    <title>coffee_genome/gene-check.py at master · drozdovapb/coffee_genome</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png">
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png">
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png">
    <meta property="fb:app_id" content="1401488693436528">

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="drozdovapb/coffee_genome" name="twitter:title" /><meta content="Contribute to coffee_genome development by creating an account on GitHub." name="twitter:description" /><meta content="https://avatars2.githubusercontent.com/u/4599655?v=2&amp;s=400" name="twitter:image:src" />
<meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars2.githubusercontent.com/u/4599655?v=2&amp;s=400" property="og:image" /><meta content="drozdovapb/coffee_genome" property="og:title" /><meta content="https://github.com/drozdovapb/coffee_genome" property="og:url" /><meta content="Contribute to coffee_genome development by creating an account on GitHub." property="og:description" />

      <meta name="browser-stats-url" content="/_stats">
    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035">
    <link rel="xhr-socket" href="/_sockets">
    <meta name="pjax-timeout" content="1000">

    <meta name="msapplication-TileImage" content="/windows-tile.png">
    <meta name="msapplication-TileColor" content="#ffffff">
    <meta name="selected-link" value="repo_source" data-pjax-transient>
      <meta name="google-analytics" content="UA-3769691-2">

    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="C390E7C6:2C59:1BA3F0C:5439201F" name="octolytics-dimension-request_id" /><meta content="4599655" name="octolytics-actor-id" /><meta content="drozdovapb" name="octolytics-actor-login" /><meta content="ac68d1f4627735b52d51f859cb09ceb787bcfedfa52488de1dec15abf6129f9f" name="octolytics-actor-hash" />
    <meta content="Rails, view, blob#show" name="analytics-event" />

    
    
    <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">


    <meta content="authenticity_token" name="csrf-param" />
<meta content="iaUTTHXpehaBRFSnsCzGrVYV9iPBpjE//U62jIQSgjgZB9+DQqZ69Rmp4BCMEATtfK4tOlZnOhHBwQUTNoT2qA==" name="csrf-token" />

    <link href="https://assets-cdn.github.com/assets/github-043670bf5d45762c99c890603216d8776470fa11262837b5ba8ca37f4175d357.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://assets-cdn.github.com/assets/github2-f97cae5c72db1b1729daa66251ec6bbfed848d4af992c2f4842aed69d5cc5277.css" media="all" rel="stylesheet" type="text/css" />
    
    


    <meta http-equiv="x-pjax-version" content="1bc26f375a269a3b5dc66c8831eeb112">

      
  <meta name="description" content="Contribute to coffee_genome development by creating an account on GitHub.">
  <meta name="go-import" content="github.com/drozdovapb/coffee_genome git https://github.com/drozdovapb/coffee_genome.git">

  <meta content="4599655" name="octolytics-dimension-user_id" /><meta content="drozdovapb" name="octolytics-dimension-user_login" /><meta content="24267260" name="octolytics-dimension-repository_id" /><meta content="drozdovapb/coffee_genome" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="24267260" name="octolytics-dimension-repository_network_root_id" /><meta content="drozdovapb/coffee_genome" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/drozdovapb/coffee_genome/commits/master.atom" rel="alternate" title="Recent Commits to coffee_genome:master" type="application/atom+xml">

  </head>


  <body class="logged_in  env-production linux vis-public page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div class="wrapper">
      
      
      
      


      <div class="header header-logged-in true" role="banner">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/" data-hotkey="g d" aria-label="Homepage" ga-data-click="Header, go to dashboard, icon:logo">
  <span class="mega-octicon octicon-mark-github"></span>
</a>


      <div class="site-search repo-scope js-site-search" role="search">
          <form accept-charset="UTF-8" action="/drozdovapb/coffee_genome/search" class="js-site-search-form" data-global-search-url="/search" data-repo-search-url="/drozdovapb/coffee_genome/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
  <input type="text"
    class="js-site-search-field is-clearable"
    data-hotkey="s"
    name="q"
    placeholder="Search"
    data-global-scope-placeholder="Search GitHub"
    data-repo-scope-placeholder="Search"
    tabindex="1"
    autocapitalize="off">
  <div class="scope-badge">This repository</div>
</form>
      </div>
      <ul class="header-nav left" role="navigation">
        <li class="header-nav-item explore">
          <a class="header-nav-link" href="/explore" data-ga-click="Header, go to explore, text:explore">Explore</a>
        </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="https://gist.github.com" data-ga-click="Header, go to gist, text:gist">Gist</a>
          </li>
          <li class="header-nav-item">
            <a class="header-nav-link" href="/blog" data-ga-click="Header, go to blog, text:blog">Blog</a>
          </li>
        <li class="header-nav-item">
          <a class="header-nav-link" href="https://help.github.com" data-ga-click="Header, go to help, text:help">Help</a>
        </li>
      </ul>

    
<ul class="header-nav user-nav right" id="user-links">
  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link name" href="/drozdovapb" data-ga-click="Header, go to profile, text:username">
      <img alt="drozdovapb" class="avatar" data-user="4599655" height="20" src="https://avatars0.githubusercontent.com/u/4599655?v=2&amp;s=40" width="20" />
      <span class="css-truncate">
        <span class="css-truncate-target">drozdovapb</span>
      </span>
    </a>
  </li>

  <li class="header-nav-item dropdown js-menu-container">
    <a class="header-nav-link js-menu-target tooltipped tooltipped-s" href="#" aria-label="Create new..." data-ga-click="Header, create new, icon:add">
      <span class="octicon octicon-plus"></span>
      <span class="dropdown-caret"></span>
    </a>

    <div class="dropdown-menu-content js-menu-content">
      
<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>


    <li class="dropdown-divider"></li>
    <li class="dropdown-header">
      <span title="drozdovapb/coffee_genome">This repository</span>
    </li>
      <li>
        <a href="/drozdovapb/coffee_genome/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
      </li>
      <li>
        <a href="/drozdovapb/coffee_genome/settings/collaboration"><span class="octicon octicon-person"></span> New collaborator</a>
      </li>
</ul>

    </div>
  </li>

  <li class="header-nav-item">
        <a href="/notifications" aria-label="You have unread notifications" class="header-nav-link notification-indicator tooltipped tooltipped-s" data-ga-click="Header, go to notifications, icon:unread" data-hotkey="g n">
        <span class="mail-status unread"></span>
        <span class="octicon octicon-inbox"></span>
</a>
  </li>

  <li class="header-nav-item">
    <a class="header-nav-link tooltipped tooltipped-s" href="/settings/profile" id="account_settings" aria-label="Settings" data-ga-click="Header, go to settings, icon:settings">
      <span class="octicon octicon-gear"></span>
    </a>
  </li>

  <li class="header-nav-item">
    <form accept-charset="UTF-8" action="/logout" class="logout-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="hM3Pn+35OttIgCukI9JgoIpzsBpWn7JmnAL0vmktDeoHoS4il7HGpkx0tiVNjKFabWcvRG7XwxTp9fM4Rf5o3g==" /></div>
      <button class="header-nav-link sign-out-button tooltipped tooltipped-s" aria-label="Sign out" data-ga-click="Header, sign out, icon:logout">
        <span class="octicon octicon-sign-out"></span>
      </button>
</form>  </li>

</ul>


    
  </div>
</div>

      

        


      <div id="start-of-content" class="accessibility-aid"></div>
          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    <div id="js-flash-container">
      
    </div>
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        
<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="qUjaY4nP34GJ7FTfP75YCeubRetcc/+DFCGRnXCJ520DUUovKxhMuQfv4TD5pp7f5gJUoTagoBC9icCTh4KBDA==" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="24267260" />

    <div class="select-menu js-menu-container js-select-menu">
      <a class="social-count js-social-count" href="/drozdovapb/coffee_genome/watchers">
        2
      </a>
      <a href="/drozdovapb/coffee_genome/subscription"
        class="minibutton select-menu-button with-count js-menu-target" role="button" tabindex="0" aria-haspopup="true">
        <span class="js-select-button">
          <span class="octicon octicon-eye"></span>
          Unwatch
        </span>
      </a>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content" aria-hidden="true">
          <div class="select-menu-header">
            <span class="select-menu-title">Notifications</span>
            <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container" role="menu">

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">Be notified when participating or @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">Be notified of all conversations.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">Never be notified.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
    
  <div class="js-toggler-container js-social-container starring-container ">

    <form accept-charset="UTF-8" action="/drozdovapb/coffee_genome/unstar" class="js-toggler-form starred js-unstar-button" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="RRZW602lTjwSfLPYnmPVYhm4QDO4z8Kdk95HaS7d0BFtNzZ/xHVfSg6SuPZobO0eCrdVEgSXO5wI9VMDgPYA1Q==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Unstar this repository" title="Unstar drozdovapb/coffee_genome">
        <span class="octicon octicon-star"></span>
        Unstar
      </button>
        <a class="social-count js-social-count" href="/drozdovapb/coffee_genome/stargazers">
          0
        </a>
</form>
    <form accept-charset="UTF-8" action="/drozdovapb/coffee_genome/star" class="js-toggler-form unstarred js-star-button" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="0ZimQb7PkqySUyqxlZ3rqGzpImTP7nPYctxMabVMbCeue7bFeyW6R+a9wZ+QsA6xFRNfHEGnyOMzWUvHxp6GcQ==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Star this repository" title="Star drozdovapb/coffee_genome">
        <span class="octicon octicon-star"></span>
        Star
      </button>
        <a class="social-count js-social-count" href="/drozdovapb/coffee_genome/stargazers">
          0
        </a>
</form>  </div>

  </li>


        <li>
          <a href="/drozdovapb/coffee_genome/fork" class="minibutton with-count js-toggler-target fork-button tooltipped-n" title="Fork your own copy of drozdovapb/coffee_genome to your account" aria-label="Fork your own copy of drozdovapb/coffee_genome to your account" rel="facebox nofollow">
            <span class="octicon octicon-repo-forked"></span>
            Fork
          </a>
          <a href="/drozdovapb/coffee_genome/network" class="social-count">0</a>
        </li>

</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="mega-octicon octicon-repo"></span>
          <span class="author"><a href="/drozdovapb" class="url fn" itemprop="url" rel="author"><span itemprop="title">drozdovapb</span></a></span><!--
       --><span class="path-divider">/</span><!--
       --><strong><a href="/drozdovapb/coffee_genome" class="js-current-repository js-repo-home-link">coffee_genome</a></strong>

          <span class="page-context-loader">
            <img alt="" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline  ">
        <div class="repository-sidebar clearfix">
            
<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders" role="navigation" data-issue-count-url="/drozdovapb/coffee_genome/issues/counts">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped tooltipped-w" aria-label="Code">
        <a href="/drozdovapb/coffee_genome" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-hotkey="g c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /drozdovapb/coffee_genome">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped tooltipped-w" aria-label="Issues">
          <a href="/drozdovapb/coffee_genome/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g i" data-selected-links="repo_issues repo_labels repo_milestones /drozdovapb/coffee_genome/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class="js-issue-replace-counter"></span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped tooltipped-w" aria-label="Pull Requests">
        <a href="/drozdovapb/coffee_genome/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g p" data-selected-links="repo_pulls /drozdovapb/coffee_genome/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class="js-pull-replace-counter"></span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped tooltipped-w" aria-label="Wiki">
          <a href="/drozdovapb/coffee_genome/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g w" data-selected-links="repo_wiki /drozdovapb/coffee_genome/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped tooltipped-w" aria-label="Pulse">
        <a href="/drozdovapb/coffee_genome/pulse/weekly" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /drozdovapb/coffee_genome/pulse/weekly">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Graphs">
        <a href="/drozdovapb/coffee_genome/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /drozdovapb/coffee_genome/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


      <div class="sunken-menu-separator"></div>
      <ul class="sunken-menu-group">
        <li class="tooltipped tooltipped-w" aria-label="Settings">
          <a href="/drozdovapb/coffee_genome/settings" aria-label="Settings" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_settings /drozdovapb/coffee_genome/settings">
            <span class="octicon octicon-tools"></span> <span class="full-word">Settings</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
      </ul>
  </div>
</div>

              <div class="only-with-full-nav">
                
  
<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><span class="text-emphasized">HTTPS</span> clone URL</h3>
  <div class="input-group">
    <input type="text" class="input-mini input-monospace js-url-field"
           value="https://github.com/drozdovapb/coffee_genome.git" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/drozdovapb/coffee_genome.git" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><span class="text-emphasized">SSH</span> clone URL</h3>
  <div class="input-group">
    <input type="text" class="input-mini input-monospace js-url-field"
           value="git@github.com:drozdovapb/coffee_genome.git" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="git@github.com:drozdovapb/coffee_genome.git" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  
<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><span class="text-emphasized">Subversion</span> checkout URL</h3>
  <div class="input-group">
    <input type="text" class="input-mini input-monospace js-url-field"
           value="https://github.com/drozdovapb/coffee_genome" readonly="readonly">
    <span class="input-group-button">
      <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/drozdovapb/coffee_genome" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <a href="https://help.github.com/articles/which-remote-url-should-i-use" class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <span class="octicon octicon-question"></span>
  </a>
</p>



                <a href="/drozdovapb/coffee_genome/archive/master.zip"
                   class="minibutton sidebar-button"
                   aria-label="Download the contents of drozdovapb/coffee_genome as a zip file"
                   title="Download the contents of drozdovapb/coffee_genome as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          

<a href="/drozdovapb/coffee_genome/blob/52ae71890d6f9e20f6303e5ec32d8e1627deed28/Gene-check%20python%203.4%20script/gene-check.py" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:408d877732973656ef7d1a492e88a6dd -->

<div class="file-navigation">
  
<div class="select-menu js-menu-container js-select-menu left">
  <span class="minibutton select-menu-button js-menu-target css-truncate" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    title="master"
    role="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button css-truncate-target">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close" role="button" aria-label="Close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/drozdovapb/coffee_genome/blob/master/Gene-check%20python%203.4%20script/gene-check.py"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/drozdovapb/coffee_genome/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="HsTawY/bguUxtl51s8aIQ8eePdsXOBpWTsLtyIrVw7UGGuQAEX4KMlA+pTdGgsIUEOmqH7GnYSU8LHDAg1gJcw==" /></div>
            <span class="octicon octicon-git-branch select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master">
            <input type="hidden" name="path" id="path" value="Gene-check python 3.4 script/gene-check.py">
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="button-group right">
    <a href="/drozdovapb/coffee_genome/find/master"
          class="js-show-file-finder minibutton empty-icon tooltipped tooltipped-s"
          data-pjax
          data-hotkey="t"
          aria-label="Quickly jump between files">
      <span class="octicon octicon-list-unordered"></span>
    </a>
    <button class="js-zeroclipboard minibutton zeroclipboard-button"
          data-clipboard-text="Gene-check python 3.4 script/gene-check.py"
          aria-label="Copy to clipboard"
          data-copied-hint="Copied!">
      <span class="octicon octicon-clippy"></span>
    </button>
  </div>

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/drozdovapb/coffee_genome" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">coffee_genome</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/drozdovapb/coffee_genome/tree/master/Gene-check%20python%203.4%20script" class="" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Gene-check python 3.4 script</span></a></span><span class="separator"> / </span><strong class="final-path">gene-check.py</strong>
  </div>
</div>


  <div class="commit file-history-tease">
    <div class="file-history-tease-header">
        <img alt="Alexandr" class="avatar" data-user="5871853" height="24" src="https://avatars3.githubusercontent.com/u/5871853?v=2&amp;s=48" width="24" />
        <span class="author"><a href="/Whiteha" rel="contributor">Whiteha</a></span>
        <time datetime="2014-10-10T23:03:52Z" is="relative-time">Oct 11, 2014</time>
        <div class="commit-title">
            <a href="/drozdovapb/coffee_genome/commit/52ae71890d6f9e20f6303e5ec32d8e1627deed28" class="message" data-pjax="true" title="Add exonerate parser to gene-check">Add exonerate parser to gene-check</a>
        </div>
    </div>

    <div class="participation">
      <p class="quickstat">
        <a href="#blob_contributors_box" rel="facebox">
          <strong>1</strong>
           contributor
        </a>
      </p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img alt="Alexandr" data-user="5871853" height="24" src="https://avatars3.githubusercontent.com/u/5871853?v=2&amp;s=48" width="24" />
            <a href="/Whiteha">Whiteha</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file-box">
  <div class="file">
    <div class="meta clearfix">
      <div class="info file-name">
          <span>228 lines (188 sloc)</span>
          <span class="meta-divider"></span>
        <span>8.831 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
          <a href="/drozdovapb/coffee_genome/raw/master/Gene-check%20python%203.4%20script/gene-check.py" class="minibutton " id="raw-url">Raw</a>
            <a href="/drozdovapb/coffee_genome/blame/master/Gene-check%20python%203.4%20script/gene-check.py" class="minibutton js-update-url-with-hash">Blame</a>
          <a href="/drozdovapb/coffee_genome/commits/master/Gene-check%20python%203.4%20script/gene-check.py" class="minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->


              <a class="octicon-button js-update-url-with-hash"
                 href="/drozdovapb/coffee_genome/edit/master/Gene-check%20python%203.4%20script/gene-check.py"
                 data-method="post" rel="nofollow" data-hotkey="e"><span class="octicon octicon-pencil"></span></a>

            <a class="octicon-button danger"
               href="/drozdovapb/coffee_genome/delete/master/Gene-check%20python%203.4%20script/gene-check.py"
               data-method="post" data-test-id="delete-blob-file" rel="nofollow">
          <span class="octicon octicon-trashcan"></span>
        </a>
      </div><!-- /.actions -->
    </div>
    
  <div class="blob-wrapper data type-python">
      <table class="highlight tab-size-8 js-file-line-container">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code js-file-line"><span class="n">__author__</span> <span class="o">=</span> <span class="s">&#39;anonym&#39;</span></td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code js-file-line"><span class="kn">import</span> <span class="nn">argparse</span></td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code js-file-line"><span class="kn">import</span> <span class="nn">io</span></td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code js-file-line"><span class="c">#################### Глобальные переменные ####################</span></td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code js-file-line"><span class="c"># DEBUG_ выставляйте в True при дебаге</span></td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code js-file-line"><span class="n">DEBUG_</span> <span class="o">=</span> <span class="bp">True</span></td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code js-file-line"><span class="n">codon_start_dna</span> <span class="o">=</span> <span class="s">&#39;ATG&#39;</span></td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code js-file-line"><span class="n">codon_stop_dna</span>  <span class="o">=</span> <span class="p">[</span><span class="s">&#39;TAA&#39;</span><span class="p">,</span><span class="s">&#39;TAG&#39;</span><span class="p">,</span><span class="s">&#39;TGA&#39;</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code js-file-line"><span class="n">codon_start_rna</span> <span class="o">=</span> <span class="s">&#39;AUG&#39;</span></td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code js-file-line"><span class="n">codon_stop_rna</span>  <span class="o">=</span> <span class="p">[</span><span class="s">&#39;UAA&#39;</span><span class="p">,</span><span class="s">&#39;UAG&#39;</span><span class="p">,</span><span class="s">&#39;UGA&#39;</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code js-file-line"><span class="n">RnaToProtein</span> <span class="o">=</span>\</td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code js-file-line">    <span class="p">{</span><span class="s">&#39;UUU&#39;</span><span class="p">:</span> <span class="s">&#39;F&#39;</span><span class="p">,</span> <span class="s">&#39;CUU&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;AUU&#39;</span><span class="p">:</span> <span class="s">&#39;I&#39;</span><span class="p">,</span> <span class="s">&#39;GUU&#39;</span><span class="p">:</span> <span class="s">&#39;V&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code js-file-line">     <span class="s">&#39;UUC&#39;</span><span class="p">:</span> <span class="s">&#39;F&#39;</span><span class="p">,</span> <span class="s">&#39;CUC&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;AUC&#39;</span><span class="p">:</span> <span class="s">&#39;I&#39;</span><span class="p">,</span> <span class="s">&#39;GUC&#39;</span><span class="p">:</span> <span class="s">&#39;V&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code js-file-line">     <span class="s">&#39;UUA&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;CUA&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;AUA&#39;</span><span class="p">:</span> <span class="s">&#39;I&#39;</span><span class="p">,</span> <span class="s">&#39;GUA&#39;</span><span class="p">:</span> <span class="s">&#39;V&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code js-file-line">     <span class="s">&#39;UUG&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;CUG&#39;</span><span class="p">:</span> <span class="s">&#39;L&#39;</span><span class="p">,</span> <span class="s">&#39;AUG&#39;</span><span class="p">:</span> <span class="s">&#39;M&#39;</span><span class="p">,</span> <span class="s">&#39;GUG&#39;</span><span class="p">:</span> <span class="s">&#39;V&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code js-file-line">     <span class="s">&#39;UCU&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;CCU&#39;</span><span class="p">:</span> <span class="s">&#39;P&#39;</span><span class="p">,</span> <span class="s">&#39;ACU&#39;</span><span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">,</span> <span class="s">&#39;GCU&#39;</span><span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code js-file-line">     <span class="s">&#39;UCC&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;CCC&#39;</span><span class="p">:</span> <span class="s">&#39;P&#39;</span><span class="p">,</span> <span class="s">&#39;ACC&#39;</span><span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">,</span> <span class="s">&#39;GCC&#39;</span><span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code js-file-line">     <span class="s">&#39;UCA&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;CCA&#39;</span><span class="p">:</span> <span class="s">&#39;P&#39;</span><span class="p">,</span> <span class="s">&#39;ACA&#39;</span><span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">,</span> <span class="s">&#39;GCA&#39;</span><span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code js-file-line">     <span class="s">&#39;UCG&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;CCG&#39;</span><span class="p">:</span> <span class="s">&#39;P&#39;</span><span class="p">,</span> <span class="s">&#39;ACG&#39;</span><span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">,</span> <span class="s">&#39;GCG&#39;</span><span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code js-file-line">     <span class="s">&#39;UAU&#39;</span><span class="p">:</span> <span class="s">&#39;Y&#39;</span><span class="p">,</span> <span class="s">&#39;CAU&#39;</span><span class="p">:</span> <span class="s">&#39;H&#39;</span><span class="p">,</span> <span class="s">&#39;AAU&#39;</span><span class="p">:</span> <span class="s">&#39;N&#39;</span><span class="p">,</span> <span class="s">&#39;GAU&#39;</span><span class="p">:</span> <span class="s">&#39;D&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code js-file-line">     <span class="s">&#39;UAC&#39;</span><span class="p">:</span> <span class="s">&#39;Y&#39;</span><span class="p">,</span> <span class="s">&#39;CAC&#39;</span><span class="p">:</span> <span class="s">&#39;H&#39;</span><span class="p">,</span> <span class="s">&#39;AAC&#39;</span><span class="p">:</span> <span class="s">&#39;N&#39;</span><span class="p">,</span> <span class="s">&#39;GAC&#39;</span><span class="p">:</span> <span class="s">&#39;D&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code js-file-line">     <span class="s">&#39;UAA&#39;</span><span class="p">:</span> <span class="s">&#39;$&#39;</span><span class="p">,</span> <span class="s">&#39;CAA&#39;</span><span class="p">:</span> <span class="s">&#39;Q&#39;</span><span class="p">,</span> <span class="s">&#39;AAA&#39;</span><span class="p">:</span> <span class="s">&#39;K&#39;</span><span class="p">,</span> <span class="s">&#39;GAA&#39;</span><span class="p">:</span> <span class="s">&#39;E&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code js-file-line">     <span class="s">&#39;UAG&#39;</span><span class="p">:</span> <span class="s">&#39;$&#39;</span><span class="p">,</span> <span class="s">&#39;CAG&#39;</span><span class="p">:</span> <span class="s">&#39;Q&#39;</span><span class="p">,</span> <span class="s">&#39;AAG&#39;</span><span class="p">:</span> <span class="s">&#39;K&#39;</span><span class="p">,</span> <span class="s">&#39;GAG&#39;</span><span class="p">:</span> <span class="s">&#39;E&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code js-file-line">     <span class="s">&#39;UGU&#39;</span><span class="p">:</span> <span class="s">&#39;C&#39;</span><span class="p">,</span> <span class="s">&#39;CGU&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;AGU&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;GGU&#39;</span><span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code js-file-line">     <span class="s">&#39;UGC&#39;</span><span class="p">:</span> <span class="s">&#39;C&#39;</span><span class="p">,</span> <span class="s">&#39;CGC&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;AGC&#39;</span><span class="p">:</span> <span class="s">&#39;S&#39;</span><span class="p">,</span> <span class="s">&#39;GGC&#39;</span><span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code js-file-line">     <span class="s">&#39;UGA&#39;</span><span class="p">:</span> <span class="s">&#39;$&#39;</span><span class="p">,</span> <span class="s">&#39;CGA&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;AGA&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;GGA&#39;</span><span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code js-file-line">     <span class="s">&#39;UGG&#39;</span><span class="p">:</span> <span class="s">&#39;W&#39;</span><span class="p">,</span> <span class="s">&#39;CGG&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;AGG&#39;</span><span class="p">:</span> <span class="s">&#39;R&#39;</span><span class="p">,</span> <span class="s">&#39;GGG&#39;</span><span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">}</span></td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code js-file-line"><span class="n">ProteinToRna</span> <span class="o">=</span> <span class="p">{</span><span class="s">&#39;A&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;GCA&#39;</span><span class="p">,</span> <span class="s">&#39;GCC&#39;</span><span class="p">,</span> <span class="s">&#39;GCU&#39;</span><span class="p">,</span> <span class="s">&#39;GCG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code js-file-line">                <span class="s">&#39;N&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;AAU&#39;</span><span class="p">,</span> <span class="s">&#39;AAC&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code js-file-line">                <span class="s">&#39;F&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UUC&#39;</span><span class="p">,</span> <span class="s">&#39;UUU&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code js-file-line">                <span class="s">&#39;Y&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UAC&#39;</span><span class="p">,</span> <span class="s">&#39;UAU&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code js-file-line">                <span class="s">&#39;$&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UGA&#39;</span><span class="p">,</span> <span class="s">&#39;UAA&#39;</span><span class="p">,</span> <span class="s">&#39;UAG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code js-file-line">                <span class="s">&#39;L&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UUG&#39;</span><span class="p">,</span> <span class="s">&#39;CUA&#39;</span><span class="p">,</span> <span class="s">&#39;CUU&#39;</span><span class="p">,</span> <span class="s">&#39;CUC&#39;</span><span class="p">,</span> <span class="s">&#39;CUG&#39;</span><span class="p">,</span> <span class="s">&#39;UUA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code js-file-line">                <span class="s">&#39;I&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;AUC&#39;</span><span class="p">,</span> <span class="s">&#39;AUU&#39;</span><span class="p">,</span> <span class="s">&#39;AUA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code js-file-line">                <span class="s">&#39;P&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;CCU&#39;</span><span class="p">,</span> <span class="s">&#39;CCG&#39;</span><span class="p">,</span> <span class="s">&#39;CCC&#39;</span><span class="p">,</span> <span class="s">&#39;CCA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code js-file-line">                <span class="s">&#39;Q&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;CAG&#39;</span><span class="p">,</span> <span class="s">&#39;CAA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code js-file-line">                <span class="s">&#39;M&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;AUG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code js-file-line">                <span class="s">&#39;S&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;AGC&#39;</span><span class="p">,</span> <span class="s">&#39;AGU&#39;</span><span class="p">,</span> <span class="s">&#39;UCG&#39;</span><span class="p">,</span> <span class="s">&#39;UCU&#39;</span><span class="p">,</span> <span class="s">&#39;UCC&#39;</span><span class="p">,</span> <span class="s">&#39;UCA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code js-file-line">                <span class="s">&#39;E&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;GAA&#39;</span><span class="p">,</span> <span class="s">&#39;GAG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code js-file-line">                <span class="s">&#39;C&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UGC&#39;</span><span class="p">,</span> <span class="s">&#39;UGU&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code js-file-line">                <span class="s">&#39;K&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;AAA&#39;</span><span class="p">,</span> <span class="s">&#39;AAG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code js-file-line">                <span class="s">&#39;R&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;CGA&#39;</span><span class="p">,</span> <span class="s">&#39;CGG&#39;</span><span class="p">,</span> <span class="s">&#39;AGG&#39;</span><span class="p">,</span> <span class="s">&#39;CGC&#39;</span><span class="p">,</span> <span class="s">&#39;AGA&#39;</span><span class="p">,</span> <span class="s">&#39;CGU&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code js-file-line">                <span class="s">&#39;W&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;UGG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code js-file-line">                <span class="s">&#39;T&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;ACU&#39;</span><span class="p">,</span> <span class="s">&#39;ACC&#39;</span><span class="p">,</span> <span class="s">&#39;ACG&#39;</span><span class="p">,</span> <span class="s">&#39;ACA&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code js-file-line">                <span class="s">&#39;V&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;GUA&#39;</span><span class="p">,</span> <span class="s">&#39;GUU&#39;</span><span class="p">,</span> <span class="s">&#39;GUC&#39;</span><span class="p">,</span> <span class="s">&#39;GUG&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code js-file-line">                <span class="s">&#39;D&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;GAU&#39;</span><span class="p">,</span> <span class="s">&#39;GAC&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code js-file-line">                <span class="s">&#39;G&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;GGU&#39;</span><span class="p">,</span> <span class="s">&#39;GGA&#39;</span><span class="p">,</span> <span class="s">&#39;GGG&#39;</span><span class="p">,</span> <span class="s">&#39;GGC&#39;</span><span class="p">],</span></td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code js-file-line">                <span class="s">&#39;H&#39;</span><span class="p">:</span> <span class="p">[</span><span class="s">&#39;CAU&#39;</span><span class="p">,</span> <span class="s">&#39;CAC&#39;</span><span class="p">]}</span></td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code js-file-line"><span class="n">DnaComplimentary</span> <span class="o">=</span> <span class="p">{</span><span class="s">&#39;A&#39;</span> <span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">,</span> <span class="s">&#39;C&#39;</span> <span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span> <span class="s">&#39;G&#39;</span> <span class="p">:</span> <span class="s">&#39;C&#39;</span><span class="p">,</span> <span class="s">&#39;T&#39;</span><span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">}</span></td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code js-file-line"><span class="n">DnaToRna</span> <span class="o">=</span> <span class="p">{</span><span class="s">&#39;A&#39;</span> <span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span> <span class="s">&#39;C&#39;</span> <span class="p">:</span> <span class="s">&#39;C&#39;</span><span class="p">,</span> <span class="s">&#39;G&#39;</span> <span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span> <span class="s">&#39;T&#39;</span><span class="p">:</span> <span class="s">&#39;U&#39;</span><span class="p">}</span></td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code js-file-line"><span class="n">RnaToDna</span> <span class="o">=</span> <span class="p">{</span><span class="s">&#39;A&#39;</span> <span class="p">:</span> <span class="s">&#39;A&#39;</span><span class="p">,</span> <span class="s">&#39;C&#39;</span> <span class="p">:</span> <span class="s">&#39;C&#39;</span><span class="p">,</span> <span class="s">&#39;G&#39;</span> <span class="p">:</span> <span class="s">&#39;G&#39;</span><span class="p">,</span> <span class="s">&#39;U&#39;</span><span class="p">:</span> <span class="s">&#39;T&#39;</span><span class="p">}</span></td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code js-file-line"><span class="k">def</span> <span class="nf">GetReversComplimentaryDna</span><span class="p">(</span><span class="n">dna</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code js-file-line">    <span class="n">result</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code js-file-line">    <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="nb">len</span><span class="p">(</span><span class="n">dna</span><span class="p">)):</span></td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code js-file-line">        <span class="n">result</span> <span class="o">+=</span> <span class="n">DnaComplimentary</span><span class="p">[</span><span class="n">dna</span><span class="p">[</span><span class="n">i</span><span class="p">]]</span></td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code js-file-line">    <span class="k">return</span> <span class="n">result</span><span class="p">[::</span><span class="o">-</span><span class="mi">1</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code js-file-line"><span class="c">####################        Функции        ####################</span></td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code js-file-line"><span class="c"># Для выода диагностических сообщейний при дебаге</span></td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code js-file-line"><span class="k">def</span> <span class="nf">printdbg</span><span class="p">(</span><span class="o">*</span><span class="n">args</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code js-file-line">    <span class="k">if</span> <span class="n">DEBUG_</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code js-file-line">        <span class="k">print</span><span class="p">(</span><span class="o">*</span><span class="n">args</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code js-file-line"><span class="c">#  Временный класс, мб вообще удалю или перепишу</span></td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code js-file-line"><span class="k">class</span> <span class="nc">SimpleReport</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code js-file-line">    <span class="n">flag</span> <span class="o">=</span> <span class="mi">0</span></td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code js-file-line">    <span class="n">description</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">__init__</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">fl</span><span class="p">,</span> <span class="n">descr</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code js-file-line">        <span class="bp">self</span><span class="o">.</span><span class="n">flag</span> <span class="o">=</span> <span class="n">fl</span></td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code js-file-line">        <span class="bp">self</span><span class="o">.</span><span class="n">description</span> <span class="o">=</span> <span class="n">descr</span></td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">__repr__</span><span class="p">(</span><span class="bp">self</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code js-file-line">        <span class="k">return</span> <span class="p">[</span><span class="bp">self</span><span class="o">.</span><span class="n">flag</span><span class="p">,</span> <span class="bp">self</span><span class="o">.</span><span class="n">description</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code js-file-line"><span class="c"># Для нашей работы нам наужно сам геном и файл выравнивания с белками</span></td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code js-file-line"><span class="k">def</span> <span class="nf">TotalCheck</span><span class="p">(</span><span class="n">_genome</span><span class="p">,</span> <span class="n">_protein_alignment</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code js-file-line">    <span class="c"># Набор функций для проверки по отдельным критериям.</span></td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code js-file-line">    <span class="c"># Тк мы работаем над небольшими геномами, то кажется можем позволить себе стратегию модульной проверки</span></td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code js-file-line">    <span class="c"># те отдельной проверки результата по каждому критерию независимо. Каждая функция проверки будет возвращать</span></td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code js-file-line">    <span class="c"># подмножество из _protein_alignment удовлетворяющее ее критерию.</span></td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code js-file-line">    <span class="c"># Итоговый результат должен быть объединением всех промежуточных проверок - это и будет список надежных генов.</span></td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">CheckOpenReadingFrame</span><span class="p">(</span><span class="n">genome</span><span class="p">,</span> <span class="n">protein_alignment</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code js-file-line">        <span class="k">return</span> <span class="n">SimpleReport</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="s">&quot;FAILURE&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">CheckExoneIntroneStructure</span><span class="p">(</span><span class="n">genome</span><span class="p">,</span> <span class="n">protein_alignment</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code js-file-line">        <span class="k">return</span> <span class="n">SimpleReport</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="s">&quot;FAILURE&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">CheckRibosomalBindingSite</span><span class="p">(</span><span class="n">genome</span><span class="p">,</span> <span class="n">protein_alignment</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code js-file-line">        <span class="k">return</span> <span class="n">SimpleReport</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="s">&quot;FAILURE&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code js-file-line">    <span class="k">def</span> <span class="nf">CheckCodonUsageBias</span><span class="p">(</span><span class="n">genome</span><span class="p">,</span> <span class="n">protein_alignment</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code js-file-line">        <span class="k">return</span> <span class="n">SimpleReport</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="s">&quot;FAILURE&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code js-file-line">    <span class="n">check_stack</span> <span class="o">=</span> <span class="p">[</span><span class="n">CheckOpenReadingFrame</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code js-file-line">                  <span class="n">CheckExoneIntroneStructure</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code js-file-line">                  <span class="n">CheckRibosomalBindingSite</span><span class="p">,</span></td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code js-file-line">                  <span class="n">CheckCodonUsageBias</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L112" class="blob-num js-line-number" data-line-number="112"></td>
        <td id="LC112" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L113" class="blob-num js-line-number" data-line-number="113"></td>
        <td id="LC113" class="blob-code js-file-line">    <span class="c"># В разработке</span></td>
      </tr>
      <tr>
        <td id="L114" class="blob-num js-line-number" data-line-number="114"></td>
        <td id="LC114" class="blob-code js-file-line">    <span class="n">report</span> <span class="o">=</span> <span class="p">[]</span></td>
      </tr>
      <tr>
        <td id="L115" class="blob-num js-line-number" data-line-number="115"></td>
        <td id="LC115" class="blob-code js-file-line">    <span class="k">for</span> <span class="n">function</span> <span class="ow">in</span> <span class="n">check_stack</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L116" class="blob-num js-line-number" data-line-number="116"></td>
        <td id="LC116" class="blob-code js-file-line">        <span class="n">report</span><span class="o">.</span><span class="n">append</span><span class="p">(</span><span class="n">function</span><span class="p">(</span><span class="n">_genome</span><span class="p">,</span> <span class="n">_protein_alignment</span><span class="p">))</span></td>
      </tr>
      <tr>
        <td id="L117" class="blob-num js-line-number" data-line-number="117"></td>
        <td id="LC117" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L118" class="blob-num js-line-number" data-line-number="118"></td>
        <td id="LC118" class="blob-code js-file-line">    <span class="k">return</span> <span class="n">report</span></td>
      </tr>
      <tr>
        <td id="L119" class="blob-num js-line-number" data-line-number="119"></td>
        <td id="LC119" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L120" class="blob-num js-line-number" data-line-number="120"></td>
        <td id="LC120" class="blob-code js-file-line"><span class="c"># {имя_белка: [(start_pos, stop_pos), ...]}</span></td>
      </tr>
      <tr>
        <td id="L121" class="blob-num js-line-number" data-line-number="121"></td>
        <td id="LC121" class="blob-code js-file-line"><span class="k">def</span> <span class="nf">ReadProteinAlignment</span><span class="p">(</span><span class="n">filename</span> <span class="o">=</span> <span class="s">&quot;&quot;</span><span class="p">,</span> <span class="n">file_from</span> <span class="o">=</span> <span class="s">&quot;exonerate&quot;</span><span class="p">):</span></td>
      </tr>
      <tr>
        <td id="L122" class="blob-num js-line-number" data-line-number="122"></td>
        <td id="LC122" class="blob-code js-file-line">    <span class="k">if</span> <span class="n">filename</span> <span class="o">==</span> <span class="s">&quot;&quot;</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L123" class="blob-num js-line-number" data-line-number="123"></td>
        <td id="LC123" class="blob-code js-file-line">        <span class="ne">Exception</span><span class="p">(</span><span class="s">&quot;Set the file name for check&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L124" class="blob-num js-line-number" data-line-number="124"></td>
        <td id="LC124" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L125" class="blob-num js-line-number" data-line-number="125"></td>
        <td id="LC125" class="blob-code js-file-line">    <span class="n">file_from</span> <span class="o">=</span> <span class="n">file_from</span><span class="o">.</span><span class="n">lower</span><span class="p">()</span></td>
      </tr>
      <tr>
        <td id="L126" class="blob-num js-line-number" data-line-number="126"></td>
        <td id="LC126" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L127" class="blob-num js-line-number" data-line-number="127"></td>
        <td id="LC127" class="blob-code js-file-line">    <span class="nb">file</span> <span class="o">=</span> <span class="nb">open</span><span class="p">(</span><span class="n">filename</span><span class="p">,</span> <span class="s">&#39;r&#39;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L128" class="blob-num js-line-number" data-line-number="128"></td>
        <td id="LC128" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L129" class="blob-num js-line-number" data-line-number="129"></td>
        <td id="LC129" class="blob-code js-file-line">    <span class="n">str_file</span> <span class="o">=</span> <span class="nb">file</span><span class="o">.</span><span class="n">read</span><span class="p">()</span></td>
      </tr>
      <tr>
        <td id="L130" class="blob-num js-line-number" data-line-number="130"></td>
        <td id="LC130" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L131" class="blob-num js-line-number" data-line-number="131"></td>
        <td id="LC131" class="blob-code js-file-line">    <span class="n">result</span> <span class="o">=</span> <span class="p">{}</span></td>
      </tr>
      <tr>
        <td id="L132" class="blob-num js-line-number" data-line-number="132"></td>
        <td id="LC132" class="blob-code js-file-line">    <span class="c"># В разработке, поэтому пишем парсер спрева только для exonerate&#39;a</span></td>
      </tr>
      <tr>
        <td id="L133" class="blob-num js-line-number" data-line-number="133"></td>
        <td id="LC133" class="blob-code js-file-line">    <span class="k">if</span> <span class="n">file_from</span> <span class="o">==</span> <span class="s">&quot;exonerate&quot;</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L134" class="blob-num js-line-number" data-line-number="134"></td>
        <td id="LC134" class="blob-code js-file-line">        <span class="c"># Парсим exonerate</span></td>
      </tr>
      <tr>
        <td id="L135" class="blob-num js-line-number" data-line-number="135"></td>
        <td id="LC135" class="blob-code js-file-line">        <span class="n">first_token</span> <span class="o">=</span> <span class="s">&quot;Query:&quot;</span></td>
      </tr>
      <tr>
        <td id="L136" class="blob-num js-line-number" data-line-number="136"></td>
        <td id="LC136" class="blob-code js-file-line">        <span class="n">second_token</span> <span class="o">=</span> <span class="s">&quot;Target range:&quot;</span></td>
      </tr>
      <tr>
        <td id="L137" class="blob-num js-line-number" data-line-number="137"></td>
        <td id="LC137" class="blob-code js-file-line">        <span class="n">offset</span> <span class="o">=</span> <span class="mi">1</span></td>
      </tr>
      <tr>
        <td id="L138" class="blob-num js-line-number" data-line-number="138"></td>
        <td id="LC138" class="blob-code js-file-line">        <span class="n">pos</span> <span class="o">=</span> <span class="n">offset</span></td>
      </tr>
      <tr>
        <td id="L139" class="blob-num js-line-number" data-line-number="139"></td>
        <td id="LC139" class="blob-code js-file-line">        <span class="n">protein_name</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L140" class="blob-num js-line-number" data-line-number="140"></td>
        <td id="LC140" class="blob-code js-file-line">        <span class="n">left_adr</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L141" class="blob-num js-line-number" data-line-number="141"></td>
        <td id="LC141" class="blob-code js-file-line">        <span class="n">right_adr</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L142" class="blob-num js-line-number" data-line-number="142"></td>
        <td id="LC142" class="blob-code js-file-line">        <span class="k">while</span> <span class="bp">True</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L143" class="blob-num js-line-number" data-line-number="143"></td>
        <td id="LC143" class="blob-code js-file-line">            <span class="c"># Получаем имя белка</span></td>
      </tr>
      <tr>
        <td id="L144" class="blob-num js-line-number" data-line-number="144"></td>
        <td id="LC144" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">=</span> <span class="n">str_file</span><span class="o">.</span><span class="n">find</span><span class="p">(</span><span class="n">first_token</span><span class="p">,</span> <span class="n">pos</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L145" class="blob-num js-line-number" data-line-number="145"></td>
        <td id="LC145" class="blob-code js-file-line">            <span class="k">if</span> <span class="n">pos</span> <span class="o">==</span> <span class="o">-</span><span class="mi">1</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L146" class="blob-num js-line-number" data-line-number="146"></td>
        <td id="LC146" class="blob-code js-file-line">                <span class="k">break</span></td>
      </tr>
      <tr>
        <td id="L147" class="blob-num js-line-number" data-line-number="147"></td>
        <td id="LC147" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">+=</span> <span class="nb">len</span><span class="p">(</span><span class="n">first_token</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L148" class="blob-num js-line-number" data-line-number="148"></td>
        <td id="LC148" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">+=</span> <span class="n">offset</span></td>
      </tr>
      <tr>
        <td id="L149" class="blob-num js-line-number" data-line-number="149"></td>
        <td id="LC149" class="blob-code js-file-line">            <span class="n">protein_name</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L150" class="blob-num js-line-number" data-line-number="150"></td>
        <td id="LC150" class="blob-code js-file-line">            <span class="k">while</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span> <span class="o">!=</span> <span class="s">&#39;</span><span class="se">\n</span><span class="s">&#39;</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L151" class="blob-num js-line-number" data-line-number="151"></td>
        <td id="LC151" class="blob-code js-file-line">                <span class="n">protein_name</span> <span class="o">+=</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L152" class="blob-num js-line-number" data-line-number="152"></td>
        <td id="LC152" class="blob-code js-file-line">                <span class="n">pos</span> <span class="o">+=</span> <span class="mi">1</span></td>
      </tr>
      <tr>
        <td id="L153" class="blob-num js-line-number" data-line-number="153"></td>
        <td id="LC153" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L154" class="blob-num js-line-number" data-line-number="154"></td>
        <td id="LC154" class="blob-code js-file-line">            <span class="c"># Получаем адреса выравнивания</span></td>
      </tr>
      <tr>
        <td id="L155" class="blob-num js-line-number" data-line-number="155"></td>
        <td id="LC155" class="blob-code js-file-line">            <span class="n">left_adr</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L156" class="blob-num js-line-number" data-line-number="156"></td>
        <td id="LC156" class="blob-code js-file-line">            <span class="n">right_adr</span> <span class="o">=</span> <span class="s">&quot;&quot;</span></td>
      </tr>
      <tr>
        <td id="L157" class="blob-num js-line-number" data-line-number="157"></td>
        <td id="LC157" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">=</span> <span class="n">str_file</span><span class="o">.</span><span class="n">find</span><span class="p">(</span><span class="n">second_token</span><span class="p">,</span> <span class="n">pos</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L158" class="blob-num js-line-number" data-line-number="158"></td>
        <td id="LC158" class="blob-code js-file-line">            <span class="k">if</span> <span class="n">pos</span> <span class="o">==</span> <span class="o">-</span><span class="mi">1</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L159" class="blob-num js-line-number" data-line-number="159"></td>
        <td id="LC159" class="blob-code js-file-line">                <span class="k">break</span></td>
      </tr>
      <tr>
        <td id="L160" class="blob-num js-line-number" data-line-number="160"></td>
        <td id="LC160" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">+=</span> <span class="nb">len</span><span class="p">(</span><span class="n">second_token</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L161" class="blob-num js-line-number" data-line-number="161"></td>
        <td id="LC161" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">+=</span> <span class="n">offset</span></td>
      </tr>
      <tr>
        <td id="L162" class="blob-num js-line-number" data-line-number="162"></td>
        <td id="LC162" class="blob-code js-file-line">            <span class="c"># Левый</span></td>
      </tr>
      <tr>
        <td id="L163" class="blob-num js-line-number" data-line-number="163"></td>
        <td id="LC163" class="blob-code js-file-line">            <span class="k">while</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span> <span class="o">!=</span> <span class="s">&#39; &#39;</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L164" class="blob-num js-line-number" data-line-number="164"></td>
        <td id="LC164" class="blob-code js-file-line">                <span class="n">left_adr</span> <span class="o">+=</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L165" class="blob-num js-line-number" data-line-number="165"></td>
        <td id="LC165" class="blob-code js-file-line">                <span class="n">pos</span> <span class="o">+=</span> <span class="mi">1</span></td>
      </tr>
      <tr>
        <td id="L166" class="blob-num js-line-number" data-line-number="166"></td>
        <td id="LC166" class="blob-code js-file-line">            <span class="n">pos</span> <span class="o">+=</span> <span class="nb">len</span><span class="p">(</span><span class="s">&quot;-&gt; &quot;</span><span class="p">)</span> <span class="o">+</span> <span class="mi">1</span></td>
      </tr>
      <tr>
        <td id="L167" class="blob-num js-line-number" data-line-number="167"></td>
        <td id="LC167" class="blob-code js-file-line">            <span class="c"># Правый</span></td>
      </tr>
      <tr>
        <td id="L168" class="blob-num js-line-number" data-line-number="168"></td>
        <td id="LC168" class="blob-code js-file-line">            <span class="k">while</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span> <span class="o">!=</span> <span class="s">&#39; &#39;</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L169" class="blob-num js-line-number" data-line-number="169"></td>
        <td id="LC169" class="blob-code js-file-line">                <span class="n">right_adr</span> <span class="o">+=</span> <span class="n">str_file</span><span class="p">[</span><span class="n">pos</span><span class="p">]</span></td>
      </tr>
      <tr>
        <td id="L170" class="blob-num js-line-number" data-line-number="170"></td>
        <td id="LC170" class="blob-code js-file-line">                <span class="n">pos</span> <span class="o">+=</span> <span class="mi">1</span></td>
      </tr>
      <tr>
        <td id="L171" class="blob-num js-line-number" data-line-number="171"></td>
        <td id="LC171" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L172" class="blob-num js-line-number" data-line-number="172"></td>
        <td id="LC172" class="blob-code js-file-line">            <span class="k">if</span> <span class="n">protein_name</span> <span class="ow">in</span> <span class="n">result</span><span class="o">.</span><span class="n">keys</span><span class="p">():</span></td>
      </tr>
      <tr>
        <td id="L173" class="blob-num js-line-number" data-line-number="173"></td>
        <td id="LC173" class="blob-code js-file-line">                <span class="n">result</span><span class="p">[</span><span class="n">protein_name</span><span class="p">]</span><span class="o">.</span><span class="n">append</span><span class="p">([</span><span class="nb">int</span><span class="p">(</span><span class="n">left_adr</span><span class="p">),</span> <span class="nb">int</span><span class="p">(</span><span class="n">right_adr</span><span class="p">)])</span></td>
      </tr>
      <tr>
        <td id="L174" class="blob-num js-line-number" data-line-number="174"></td>
        <td id="LC174" class="blob-code js-file-line">            <span class="k">else</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L175" class="blob-num js-line-number" data-line-number="175"></td>
        <td id="LC175" class="blob-code js-file-line">                <span class="n">result</span><span class="p">[</span><span class="n">protein_name</span><span class="p">]</span> <span class="o">=</span> <span class="p">[</span><span class="nb">int</span><span class="p">(</span><span class="n">left_adr</span><span class="p">),</span> <span class="nb">int</span><span class="p">(</span><span class="n">right_adr</span><span class="p">)]</span></td>
      </tr>
      <tr>
        <td id="L176" class="blob-num js-line-number" data-line-number="176"></td>
        <td id="LC176" class="blob-code js-file-line">    <span class="k">else</span><span class="p">:</span></td>
      </tr>
      <tr>
        <td id="L177" class="blob-num js-line-number" data-line-number="177"></td>
        <td id="LC177" class="blob-code js-file-line">        <span class="ne">Exception</span><span class="p">(</span><span class="s">&quot;Now available only for the exonerate&quot;</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L178" class="blob-num js-line-number" data-line-number="178"></td>
        <td id="LC178" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L179" class="blob-num js-line-number" data-line-number="179"></td>
        <td id="LC179" class="blob-code js-file-line">    <span class="k">return</span> <span class="n">result</span></td>
      </tr>
      <tr>
        <td id="L180" class="blob-num js-line-number" data-line-number="180"></td>
        <td id="LC180" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L181" class="blob-num js-line-number" data-line-number="181"></td>
        <td id="LC181" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L182" class="blob-num js-line-number" data-line-number="182"></td>
        <td id="LC182" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L183" class="blob-num js-line-number" data-line-number="183"></td>
        <td id="LC183" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L184" class="blob-num js-line-number" data-line-number="184"></td>
        <td id="LC184" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L185" class="blob-num js-line-number" data-line-number="185"></td>
        <td id="LC185" class="blob-code js-file-line"><span class="c">####################    Основная логика    ####################</span></td>
      </tr>
      <tr>
        <td id="L186" class="blob-num js-line-number" data-line-number="186"></td>
        <td id="LC186" class="blob-code js-file-line"><span class="k">def</span> <span class="nf">__main__</span><span class="p">():</span></td>
      </tr>
      <tr>
        <td id="L187" class="blob-num js-line-number" data-line-number="187"></td>
        <td id="LC187" class="blob-code js-file-line">    <span class="c"># Парсим аргументы командной строки и вообще ведем себя как порядочная утилита</span></td>
      </tr>
      <tr>
        <td id="L188" class="blob-num js-line-number" data-line-number="188"></td>
        <td id="LC188" class="blob-code js-file-line">    <span class="c"># parser = argparse.ArgumentParser()</span></td>
      </tr>
      <tr>
        <td id="L189" class="blob-num js-line-number" data-line-number="189"></td>
        <td id="LC189" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L190" class="blob-num js-line-number" data-line-number="190"></td>
        <td id="LC190" class="blob-code js-file-line">    <span class="c"># parser.add_argument(&#39;-sf&#39;, required = True,  help = &#39;Path to source fasta file&#39;)</span></td>
      </tr>
      <tr>
        <td id="L191" class="blob-num js-line-number" data-line-number="191"></td>
        <td id="LC191" class="blob-code js-file-line">    <span class="c"># parser.add_argument(&#39;-a&#39;, required = True, help = &#39;Path to alignment file&#39;)</span></td>
      </tr>
      <tr>
        <td id="L192" class="blob-num js-line-number" data-line-number="192"></td>
        <td id="LC192" class="blob-code js-file-line">    <span class="c"># parser.add_argument(&#39;-pn&#39;, required = True, choices = [&quot;exonerate&quot;] , help = &#39;Program name what make alignment&#39;)</span></td>
      </tr>
      <tr>
        <td id="L193" class="blob-num js-line-number" data-line-number="193"></td>
        <td id="LC193" class="blob-code js-file-line">    <span class="c"># parser.add_argument(&#39;-d&#39;, required = False, help = &#39;Destination folder for output&#39;)</span></td>
      </tr>
      <tr>
        <td id="L194" class="blob-num js-line-number" data-line-number="194"></td>
        <td id="LC194" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L195" class="blob-num js-line-number" data-line-number="195"></td>
        <td id="LC195" class="blob-code js-file-line">    <span class="c"># parse_result = vars(parser.parse_args())</span></td>
      </tr>
      <tr>
        <td id="L196" class="blob-num js-line-number" data-line-number="196"></td>
        <td id="LC196" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L197" class="blob-num js-line-number" data-line-number="197"></td>
        <td id="LC197" class="blob-code js-file-line">    <span class="c"># printdbg(&quot;Your input:&quot;, parse_result)</span></td>
      </tr>
      <tr>
        <td id="L198" class="blob-num js-line-number" data-line-number="198"></td>
        <td id="LC198" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L199" class="blob-num js-line-number" data-line-number="199"></td>
        <td id="LC199" class="blob-code js-file-line">    <span class="c"># source_fasta   = parse_result[&quot;sf&quot;]</span></td>
      </tr>
      <tr>
        <td id="L200" class="blob-num js-line-number" data-line-number="200"></td>
        <td id="LC200" class="blob-code js-file-line">    <span class="c"># alignment_file = parse_result[&quot;a&quot;]</span></td>
      </tr>
      <tr>
        <td id="L201" class="blob-num js-line-number" data-line-number="201"></td>
        <td id="LC201" class="blob-code js-file-line">    <span class="c"># program_name   = parse_result[&quot;pn&quot;]</span></td>
      </tr>
      <tr>
        <td id="L202" class="blob-num js-line-number" data-line-number="202"></td>
        <td id="LC202" class="blob-code js-file-line">    <span class="c"># destination    = parse_result[&quot;d&quot;] if parse_result[&quot;d&quot;] else &quot;output&quot;</span></td>
      </tr>
      <tr>
        <td id="L203" class="blob-num js-line-number" data-line-number="203"></td>
        <td id="LC203" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L204" class="blob-num js-line-number" data-line-number="204"></td>
        <td id="LC204" class="blob-code js-file-line">    <span class="c"># printdbg(&quot;Real destination for output:&quot;, destination)</span></td>
      </tr>
      <tr>
        <td id="L205" class="blob-num js-line-number" data-line-number="205"></td>
        <td id="LC205" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L206" class="blob-num js-line-number" data-line-number="206"></td>
        <td id="LC206" class="blob-code js-file-line">    <span class="c"># # Общий формат данных которые мы должны получать из любого выравнивания</span></td>
      </tr>
      <tr>
        <td id="L207" class="blob-num js-line-number" data-line-number="207"></td>
        <td id="LC207" class="blob-code js-file-line">    <span class="c"># # имя белка, позиции начала и конца белка в геноме, позиции интронов</span></td>
      </tr>
      <tr>
        <td id="L208" class="blob-num js-line-number" data-line-number="208"></td>
        <td id="LC208" class="blob-code js-file-line">    <span class="c"># common_protein_alignment_data = ReadProteinAlignment(alignment_file, program_name)</span></td>
      </tr>
      <tr>
        <td id="L209" class="blob-num js-line-number" data-line-number="209"></td>
        <td id="LC209" class="blob-code js-file-line">    <span class="c"># result = TotalCheck(source_fasta, common_protein_alignment_data)</span></td>
      </tr>
      <tr>
        <td id="L210" class="blob-num js-line-number" data-line-number="210"></td>
        <td id="LC210" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L211" class="blob-num js-line-number" data-line-number="211"></td>
        <td id="LC211" class="blob-code js-file-line">    <span class="c"># try:</span></td>
      </tr>
      <tr>
        <td id="L212" class="blob-num js-line-number" data-line-number="212"></td>
        <td id="LC212" class="blob-code js-file-line">    <span class="c"># #   Сохраняем отчет о проверке</span></td>
      </tr>
      <tr>
        <td id="L213" class="blob-num js-line-number" data-line-number="213"></td>
        <td id="LC213" class="blob-code js-file-line">    <span class="c">#     output_file = open(destination, &quot;w&quot;)</span></td>
      </tr>
      <tr>
        <td id="L214" class="blob-num js-line-number" data-line-number="214"></td>
        <td id="LC214" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L215" class="blob-num js-line-number" data-line-number="215"></td>
        <td id="LC215" class="blob-code js-file-line">    <span class="c">#     printdbg(result)</span></td>
      </tr>
      <tr>
        <td id="L216" class="blob-num js-line-number" data-line-number="216"></td>
        <td id="LC216" class="blob-code js-file-line">    <span class="c">#</span></td>
      </tr>
      <tr>
        <td id="L217" class="blob-num js-line-number" data-line-number="217"></td>
        <td id="LC217" class="blob-code js-file-line">    <span class="c">#     output_file.close()</span></td>
      </tr>
      <tr>
        <td id="L218" class="blob-num js-line-number" data-line-number="218"></td>
        <td id="LC218" class="blob-code js-file-line">    <span class="c"># except:</span></td>
      </tr>
      <tr>
        <td id="L219" class="blob-num js-line-number" data-line-number="219"></td>
        <td id="LC219" class="blob-code js-file-line">    <span class="c">#     print(&quot;Can&#39;t save a result to &quot; + destination)</span></td>
      </tr>
      <tr>
        <td id="L220" class="blob-num js-line-number" data-line-number="220"></td>
        <td id="LC220" class="blob-code js-file-line">    <span class="n">alignment_file</span> <span class="o">=</span> <span class="s">&quot;exonerateout&quot;</span></td>
      </tr>
      <tr>
        <td id="L221" class="blob-num js-line-number" data-line-number="221"></td>
        <td id="LC221" class="blob-code js-file-line">    <span class="n">program_name</span> <span class="o">=</span> <span class="s">&quot;exonerate&quot;</span></td>
      </tr>
      <tr>
        <td id="L222" class="blob-num js-line-number" data-line-number="222"></td>
        <td id="LC222" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L223" class="blob-num js-line-number" data-line-number="223"></td>
        <td id="LC223" class="blob-code js-file-line">    <span class="n">common_protein_alignment_data</span> <span class="o">=</span> <span class="n">ReadProteinAlignment</span><span class="p">(</span><span class="n">alignment_file</span><span class="p">,</span> <span class="n">program_name</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L224" class="blob-num js-line-number" data-line-number="224"></td>
        <td id="LC224" class="blob-code js-file-line">    <span class="k">print</span><span class="p">(</span><span class="n">common_protein_alignment_data</span><span class="p">)</span></td>
      </tr>
      <tr>
        <td id="L225" class="blob-num js-line-number" data-line-number="225"></td>
        <td id="LC225" class="blob-code js-file-line">
</td>
      </tr>
      <tr>
        <td id="L226" class="blob-num js-line-number" data-line-number="226"></td>
        <td id="LC226" class="blob-code js-file-line"><span class="c"># Запуск</span></td>
      </tr>
      <tr>
        <td id="L227" class="blob-num js-line-number" data-line-number="227"></td>
        <td id="LC227" class="blob-code js-file-line"><span class="n">__main__</span><span class="p">()</span></td>
      </tr>
</table>

  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer" role="contentinfo">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/" aria-label="Homepage">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2014 <span title="0.04740s from github-fe119-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-suggester-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="fullscreen-contents js-fullscreen-contents js-suggester-field" placeholder=""></textarea>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped tooltipped-w" aria-label="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped tooltipped-w"
      aria-label="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-x flash-close js-ajax-error-dismiss" aria-label="Dismiss error"></a>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-f4749195ce218608caf72b3ddefff5f580445386f2529c60e027cd18d1db0cb5.js" type="text/javascript"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-56901a1e60b184d90b134d3f30a8700dee7e5d313a3a70e28f6adb239d7d8797.js" type="text/javascript"></script>
      
      
        <script async src="https://www.google-analytics.com/analytics.js"></script>
  </body>
</html>

